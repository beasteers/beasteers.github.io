function softmax(x, max=0, degree=1){
  return Math.log((Math.exp(x*degree) + Math.exp(max*degree)) / 2) / degree;
}

function softmin(x, min, degree=1){
  return -softmax(-x, -min, degree);
}

function pulse(amplitude, x, mu, sigma){
  return amplitude * Math.exp(-Math.pow((x - mu) / sigma, 2))
}

var wave_equation = function(params) {
  var o = {
    beta: 1, // damping constant ish
    dt: 0.001,//0.0002, // timestep
    T0: 20, // free tension
    rho: 3, // mass per unit length
    N: 40, // number of elements to model with
    L: 1, // length of string
    FA: 500, // amplitude of max force
//     E: 0.1, // modulus of elasticity
//     l0: 0.1, // static deflection
    // r_2: 0.07111111111111111,

    
    k: 0, // the index (0|1) for the current timestep in u_k in U
    state: function(){ return this.U[o.k]; },
    prev_state: function(){ return this.U[+!o.k]; },
    
    init: function(){
      var o = this;
      this.U[!o.k] = o.fx0; // set u_k-1
      
      // initialize u_k
      for(var n = 1; n < o.N - 1; n++){
        this.U[o.k][n] = 1/2 * (
          2*(1 - o.r_2) * o.fx0[n] + o.r_2 * (o.fx0[n+1] + o.fx0[n-1]) - 
         (1 - 2*o.beta*o.dt) * o.dt * o.gx0[n]
        ) / (1 + 2*o.beta*o.dt);
      }
      return this;
    },
    
    update: function(data){
      Object.assign(this, data);
      // the distance between consecutive pts
      o.dx = o.L / o.N;
      // courant number
      o.r_2 = Math.min(0.9999, Math.pow(o.T0 / o.rho * o.dt / o.dx, 2));
      return this;
    },
    
    next: function(){
      var o = this;
      var uk = o.U[o.k]; // the current state
      var ukm1 = o.U[+!o.k]; // previous time step
      var force = this.net_force(); // the forcing function
      
      this.energy = 0;
      for(var n = 1; n < o.N - 1; n++){
        // update u_k+1
        this.U[+!o.k][n] = ( 
          2*(1 - o.r_2) * uk[n] + o.r_2 * (uk[n+1] + uk[n-1]) - // change due to the curvature of the string
         (1 - 2*o.beta*o.dt) * ukm1[n] + // relation to the last time step
         o.FA*(force[n]) * Math.pow(o.dt, 2) // force added to the system
        ) / (1 + 2*o.beta*o.dt); // decay

        this.energy += Math.pow(this.U[+!o.k][n], 2) + Math.pow((this.U[o.k][n] - this.U[+!o.k][n]) / o.dt, 2) / 2;
      }
      if(Math.log10(this.energy) < 3)
        this.stationary().stop();

      this.k  = +!o.k; // swaps the index between the two arrays
      // return u_k+1
      return this.U[+!o.k];
    },

    forces: {},
    _net_force: null,
    apply_force: function(id, force) {
      this.forces[id] = force;
      if(force && !this.interval)
        this.run();
      return this;
    },

    net_force: function(){
      return Object.values(this.forces).reduce(function(sum, force){
        if(force)
          force.forEach(function(f, i){
            if(f)
              sum[i] += f;
          });
        return sum;
      }, this._net_force.fill(0));
    },
    
    run: function(callback){
      // run the wave function at interval dt, passing the calculated state to a callback function
      var o = this; // for inside run function
      if(callback){
        this.run_func = function(){ callback(o.next()); };
      }
      else if(!this.run_func){
        this.run_func = function(){ console.log(o.next()); };
      }
      this.interval = setInterval(this.run_func, o.dt * 1000);
      this.run_func();
      return this;
    },

    stop: function(){
      // stops the equation loop
      clearInterval(this.interval);
      this.interval = null;
      return this;
    },

    stationary: function(){
      this.U[o.k].fill(0);
      return this;
    }
  };
  
  // initial conditions fx is initial position, gx is initial velocity
  o.fx0 = o.fx0 || new Array(o.N).fill(0).map(function(d, i){
    // default to a gaussian * 30cpL sin wave
    return pulse(10, i, 10, 20) * Math.sin(30*2*Math.PI * i / o.N);
  });
  o.gx0 = o.gx0 || new Array(o.N).fill(0);
  
  // storage for u_k and u_k-1
  o.U = [
    new Array(o.N).fill(0),
    new Array(o.N).fill(0)
  ];

  o._net_force = new Array(o.N).fill(0);
  
  return o.update(params || {}).init();
}


function create_wave(selector, o, params){
  // default parameters
  o = Object.assign({
    width: 900,
    height: 500,
    amplitude: 500,
    stroke_width: 3
  }, o);

  // generate wave equation
  var wave = wave_equation(params);
  console.log(wave);

  // draw canvas
  var container = d3.select(selector);
  var svg = container.append('svg')
    .attr('width', o.width)
    .attr('height', o.height)
    .attr("viewBox", "0 0 "+o.width + ' ' + o.height)
    .attr("preserveAspectRatio", "xMidYMid meet")
    .style('max-width', '100%')

  // create the string path
  var string = svg.append('path')
    .attr('stroke', 'white')
    .attr('stroke-width', o.stroke_width)
    .attr('stroke-linecap', 'round')
    .attr('fill', 'none');

  // maps value to pixel
  var x = d3.scaleLinear().domain([0, wave.N-1]).range([0, o.width]);
  var y = d3.scaleLinear()
          .clamp(true)
          .domain([-350, 350])
          .range([(o.height - o.amplitude)/2, o.height - (o.height - o.amplitude)/2]);
  
  // function to create line from data
  var line = d3.line()
    .x(function(d,i) {return x(i);})
    .y(function(d) {return y(d);})
    .curve(d3.curveNatural);


  // calculating force added to system based on mouse click position
  svg.on('mousedown', function(){
    var pos = d3.mouse(svg.node());
    var horizontal = y(0);
    var force = wave.state().map(function(pt, i){
      return pulse(2*(pos[1] - horizontal) * o.amplitude / o.height, x(i), pos[0], o.width / 20);
    });
    wave.apply_force('mouse', force);
  })
  svg.on('mouseup', function(){
    // remove the force applied after 100ms. this makes clicks slightly longer
    setTimeout(function(){
      wave.apply_force('mouse', null);
    }, 100)
  });



  var t = d3.select('#banner header').append('p')

  wave.run(function(state){
    // update line
    // var start = window.performance.now();
    string.attr('d', line(state));
    // t.text((window.performance.now() - start).toFixed(6));
    t.text(Math.log10(wave.energy));
  });

  return wave;
}


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
    beta: 5, // damping constant ish
    dt: 0.0002, // timestep
    T0: 10, // free tension
    rho: 0.5, // mass per unit length
    N: 40, // number of elements to model with
    L: 1, // length of string
    FA: 10000, // amplitude of max force
//     E: 0.1, // modulus of elasticity
//     l0: 0.1, // static deflection
    
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
    
//     frequency: function(o){
//       return ;
//     },
    
    update: function(data){
      Object.assign(this, data);
      // the distance between consecutive pts
      o.dx = o.L / o.N;
      // courant number
      o.r_2 = Math.pow(o.T0 / o.rho * o.dt / o.dx, 2);
      
      return this;
    },
    
    next: function(){
      var o = this;
      var uk = o.U[o.k]; // the current state
      var ukm1 = o.U[+!o.k]; // previous time step
      var force = o.force ? o.force() : null; // the forcing function
      
      for(var n = 1; n < o.N - 1; n++){
        // update u_k+1
        this.U[+!o.k][n] = ( 
          2*(1 - o.r_2) * uk[n] + o.r_2 * (uk[n+1] + uk[n-1]) - 
         (1 - 2*o.beta*o.dt) * ukm1[n] + o.FA*(force ? force[n] : 0) * Math.pow(o.dt, 2)
        ) / (1 + 2*o.beta*o.dt);
      }
      this.k  = +!o.k; // swaps the index between the two arrays
      // return u_k+1
      return this.U[+!o.k];
    },
    
    run: function(callback){
      var o = this; // for inside function
      var func = function(){ callback(o.next()); };
      this.interval = setInterval(func, o.dt * 1000);
      func();
    },
    
    stop: function(){
      clearInterval(this.interval);
    }
  };
  
  // initial conditions
  o.fx0 = o.fx0 || new Array(o.N).fill(0).map(function(d, i){
    return pulse(10, i, 10, 20) * Math.sin(30*2*Math.PI * i / o.N);
  });
  o.gx0 = o.gx0 || new Array(o.N).fill(0);
  o.U = [
    new Array(o.N).fill(0),
    new Array(o.N).fill(0)
  ];
  
  return o.update(params || {}).init();
}


function create_wave(selector, o, params){
  o = Object.assign({
    width: 900,
    height: 500,
    amplitude: 400,
    stroke_width: 3
  }, o);

  var wave = wave_equation(params);
  console.log(wave);

  var svg = d3.select(selector).append('svg')
  //   .attr('width', width)
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
          .domain([-200, 200])
          .range([(o.height - o.amplitude)/2, o.height - (o.height - o.amplitude)/2]);
  // function to create line from data
  var line = d3.line()
    .x(function(d,i) {return x(i);})
    .y(function(d) {return y(d);})
    .curve(d3.curveNatural);

  // var t = svg.append('text').attr('x', 10).attr('y', 0).text('hi').style('font-size', 30)


  var force = new Array(wave.N).fill(0);

  svg.on('mousedown', function(){
    var pos = d3.mouse(svg.node());
    var horizontal = y(0);
    force = wave.state().map(function(pt, i){
      var amp = 2*(pos[1] - horizontal) * o.amplitude / o.height;
      return pulse(amp, x(i), pos[0], o.width / 20);
      // return 2*(pos[1] - horizontal) * o.amplitude / o.height * Math.exp(-Math.pow((x(i) - pos[0]) / o.width * 20, 2));
    });
  })
  svg.on('mouseup', function(){
    setTimeout(function(){
      force = null;
    }, 100)
  });




  wave.update({
    force: function(){
      return force;
    }
  })

  wave.run(function(state){
    // update line
    string.attr('d', line(state));
  });

  return wave;
}



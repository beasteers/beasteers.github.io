from hyde.plugin import Plugin
from jinja2 import environmentfilter, contextfilter

# import re
import hashlib
# from slugify import slugify


def do_hash(text):
	return hashlib.md5(text.encode('utf-8')).hexdigest()

def is_list(var):
    return isinstance(var, list)


class CustomizeJinja2Plugin(Plugin):
    '''
    The curstom-filter plugin allows any
    filters added to the "filters" dictionary
    to be added to hyde
    '''

    def template_loaded(self, template):
        super(CustomizeJinja2Plugin, self).template_loaded(template)
        template.env.filters.update(dict(
        	# slugify=slugify,
        	# re=re,
        	hash=do_hash,
            is_list=is_list,
        ))
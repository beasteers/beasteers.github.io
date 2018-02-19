from hyde.ext.plugins.css import SassPlugin

import os
import glob

from fswrap import File

#
# Sass CSS
#


class SassPlugin(SassPlugin):
    _files = None
    def _should_parse_resource(self, resource):
        """
        Check user defined
        """
        content_root = self.site.config.content_root_path
        self._files = self._files or [
            os.path.relpath(f, content_root.path) # convert back to relative path
            for fp in self.settings.get('files', []) # loop thru list of files/patterns
            for f in glob.glob(content_root.child(fp)) # flatten files from file patterns
        ]
        return (resource.source_file.kind == 'scss' and 
                self.settings.get('parse', True) and 
               (not self._files or resource.relative_path in self._files))

        
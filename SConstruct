# -*- python -*-

import atexit
import os

AddOption('--production',
          dest='production',
          action='store_true',
          help='Build production versions')

AddOption('--build-pdf',
          dest='build_pdf',
          action='store_true',
          help='Build PDF file too')

AddOption('--variant',
          dest='variant',
          type='string', nargs=1, action='store',
          default='default',
          help='Build variant of the textbook')

AddOption('--chunk-size',
          dest='chunksize',
          type='int', nargs=1, action='store',
          default=250,
          help='Process this many files at once in pretex')

env = Environment(tools=['default', 'textfile', TOOL_ADD_CAT],
                  PRODUCTION=GetOption('production'),
                  VARIANT=GetOption('variant'),
                  CHUNKSIZE=GetOption('chunksize'),
                  BUILD_PDF=GetOption('build_pdf'))

env['BASE_DIR'] = env.Entry('#').get_abspath()
env['BUILD_DIR'] = env.Entry('#build').get_abspath()
env['STAGING_DIR'] = env.Entry('#staging').get_abspath()

env['UGLIFYJS'] = env.File('#/node_modules/uglify-js/bin/uglifyjs').get_abspath()
env['CLEANCSS'] = env.File('#/node_modules/clean-css-cli/bin/cleancss').get_abspath()
env['PRETEX']   = env.File('#/pretex/pretex.py').get_abspath()

build_dir   = lambda sd='': os.path.join(env['BUILD_DIR'],   sd)
staging_dir = lambda sd='': os.path.join(env['STAGING_DIR'], sd)

subpackages = []
for dirname in ['mathbook', 'mathbook-assets', 'mathbox']:
    target = dirname + '.build'
    env.Pseudo(target)
    node = env.Command(target, [], "cd " + dirname + " && scons")
    subpackages.append(node)

env.Alias('subpackages', subpackages)

env.Default('build')
env.Depends('build', 'delete_build')
env.Pseudo('delete_build')
env.Command('delete_build', [], 'rm -rf build')

if not os.path.isdir('node_modules'):
    env.Depends('build', 'node_modules')
    env.Command('node_modules', [], 'npm install')

# Create bundles

to_minify = []

ila_js = env.CatJS('$STAGING_DIR/js/ila',
                   Split('''
                   vendor/jquery.min.js
                   vendor/jquery.sticky.js
                   vendor/knowl
                   static/js/Mathbook
                   '''))
ila_css = env.CatCSS('$STAGING_DIR/css/ila',
                     Split('''
                     mathbook-assets/build/mathbook-gt
                     mathbook/css/mathbook-add-on
                     static/css/ila-add-on
                     vendor/knowlstyle
                     '''))
to_minify.append(ila_js)
to_minify.append(ila_css)

root = env.Dir('$STAGING_DIR')
for node in Flatten(to_minify):
    env.Minify('$BUILD_DIR/' + node.get_path(root), node)

# Copy static files

fonts = \
    env.Glob('mathbook-assets/stylesheets/fonts/ionicons/fonts/*') + \
    env.Glob('static/fonts/*')
for font in fonts:
    env.Command('$BUILD_DIR/css/fonts/' + os.path.basename(str(font)), font,
                Copy('$TARGET', '$SOURCE'))

for fname in ['images', 'manifest.json', 'google9ccfcae89045309c.html']:
    env.Command('$BUILD_DIR/' + fname, 'static/' + fname,
                Copy('$TARGET', '$SOURCE'))

env.SConscriptChdir(1)
env.SConscript('demos/SConscript', exports='env build_dir staging_dir')
env.SConscript('src/SConscript', exports='env build_dir staging_dir')

def when_done():
    from SCons.Script import GetBuildFailures
    if not list(GetBuildFailures()):
        print("")
        print("Build successful!  Open or reload")
        print("     http://localhost:8081/")
        print("in your browser to see the result.")

atexit.register(when_done)

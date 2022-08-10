# -*- python -*-

import atexit
import os

AddOption('--build-pdf',
          dest='build_pdf',
          action='store_true',
          help='Build PDF file too')

AddOption('--minify',
          dest='minify',
          action='store_true',
          help='Minify css and js files')

AddOption('--scratch',
          dest='scratch',
          action='store_true',
          help='Delete the build directory first')

AddOption('--delete-cache',
          dest='delete_cache',
          action='store_true',
          help='Delete all cache files.  Implies --scratch.')

AddOption('--theme',
          dest='theme',
          type='string', nargs=1, action='store',
          default='duke',
          help='Choose the theme to use (gt or duke)')

AddOption('--chunk-size',
          dest='chunksize',
          type='int', nargs=1, action='store',
          default=250,
          help='Process this many files at once in pretex')

AddOption('--production',
          dest='production',
          action='store_true',
          help='Synonym for --build-pdf --minify --scratch')

environ = dict(os.environ)
# This somehow gets set by the nix flake.  It causes <today /> to be the UNIX
# epoch.
if 'SOURCE_DATE_EPOCH' in environ:
    del environ['SOURCE_DATE_EPOCH']

env = Environment(tools=['default', 'textfile', TOOL_ADD_CAT],
                  ENV=environ,
                  BUILD_PDF=GetOption('build_pdf'),
                  MINIFY=GetOption('minify'),
                  SCRATCH=GetOption('scratch'),
                  THEME=GetOption('theme'),
                  CHUNKSIZE=GetOption('chunksize'))
if GetOption('production'):
    env['BUILD_PDF'] = True
    env['MINIFY'] = True
    env['SCRATCH'] = True

env['BASE_DIR'] = env.Entry('#').get_abspath()
env['BUILD_DIR'] = os.path.join(env['BASE_DIR'], 'build')
env['CACHE_DIR'] = os.path.join(env['BASE_DIR'], 'cache')

if GetOption('delete_cache'):
    env['SCRATCH'] = True
    env.Execute('rm -rf $CACHE_DIR')
if env['SCRATCH']:
    env.Execute('rm -rf $BUILD_DIR')

env['UGLIFYJS'] = env.File('#/node_modules/uglify-js/bin/uglifyjs').get_abspath()
env['CLEANCSS'] = env.File('#/node_modules/clean-css-cli/bin/cleancss').get_abspath()
env['PRETEX']   = env.File('#/pretex/pretex.py').get_abspath()

build_dir = lambda sd='': os.path.join(env['BUILD_DIR'], sd)
cache_dir = lambda sd='': os.path.join(env['CACHE_DIR'], sd)

subpackages = []
for dirname in ['mathbook', 'mathbook-assets', 'mathbox']:
    target = dirname + '.build'
    env.Pseudo(target)
    node = env.Command(target, [], "cd " + dirname + " && scons")
    subpackages.append(node)

env.Alias('subpackages', subpackages)

env.Alias('build-all', '$BUILD_DIR/index.html')
env.Default('build-all')

if not os.path.isdir('node_modules'):
    env.Depends('build-all', 'node_modules')
    env.Command('node_modules', [], 'npm install')

# Create bundles

to_minify = []

ila_js = env.CatJS('$CACHE_DIR/js/ila',
                   Split('''
                   vendor/jquery.min.js
                   vendor/jquery.sticky.js
                   vendor/knowl
                   static/js/Mathbook
                   '''))
ila_css = env.CatCSS('$CACHE_DIR/css/ila',
                     Split('''
                     mathbook-assets/build/mathbook-${THEME}
                     mathbook/css/mathbook-add-on
                     static/css/ila-add-on
                     static/css/ila-add-on-${THEME}
                     vendor/knowlstyle
                     '''))
to_minify.append(ila_js)
to_minify.append(ila_css)

root = env.Dir('$CACHE_DIR')
for node in Flatten(to_minify):
    path = '$BUILD_DIR/' + node.get_path(root)
    env.Minify(path, node)
    env.Depends('build-all', path)

# Copy static files

fonts = \
    env.Glob('mathbook-assets/stylesheets/fonts/ionicons/fonts/*') + \
    env.Glob('static/fonts/*')
for font in fonts:
    dep = env.Command('$BUILD_DIR/css/fonts/' + os.path.basename(str(font)), font,
                      Copy('$TARGET', '$SOURCE'))
    env.Depends('build-all', dep)

for fname in ['manifest.json', 'google9ccfcae89045309c.html']:
    dep = env.Command('$BUILD_DIR/' + fname, 'static/' + fname,
                      Copy('$TARGET', '$SOURCE'))
    env.Depends('build-all', dep)

images = \
    env.Command('$BUILD_DIR/images', 'static/images',
                Copy('$TARGET', '$SOURCE'))
env.AddPostAction(images, 'cp $BASE_DIR/static/theme-$THEME/* $BUILD_DIR/images')
env.Depends('build-all', images)

favicon = env.Command('$BUILD_DIR/favicon.ico', '$BUILD_DIR/images/logo.png',
                      Copy('$TARGET', '$SOURCE'))
env.Depends('build-all', favicon)


env.SConscriptChdir(1)
env.SConscript('demos/SConscript', exports='env build_dir cache_dir')
env.SConscript('src/SConscript', exports='env build_dir cache_dir')


def when_done():
    from SCons.Script import GetBuildFailures
    if not list(GetBuildFailures()):
        print("")
        print("Build successful!")

atexit.register(when_done)


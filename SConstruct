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

env = Environment(tools=['default', 'textfile', TOOL_ADD_CAT],
                  ENV=os.environ,
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

env.Alias('build', env['BUILD_DIR'])
env.Default('build')

if not os.path.isdir('node_modules'):
    env.Depends('build', 'node_modules')
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
    env.Minify('$BUILD_DIR/' + node.get_path(root), node)

# Copy static files

fonts = \
    env.Glob('mathbook-assets/stylesheets/fonts/ionicons/fonts/*') + \
    env.Glob('static/fonts/*')
for font in fonts:
    env.Command('$BUILD_DIR/css/fonts/' + os.path.basename(str(font)), font,
                Copy('$TARGET', '$SOURCE'))

for fname in ['manifest.json', 'google9ccfcae89045309c.html']:
    env.Command('$BUILD_DIR/' + fname, 'static/' + fname,
                Copy('$TARGET', '$SOURCE'))

images = \
    env.Command('$BUILD_DIR/images', 'static/images',
                Copy('$TARGET', '$SOURCE'))
env.AddPostAction(images, 'cp $BASE_DIR/static/theme-$THEME/* $BUILD_DIR/images')


env.SConscriptChdir(1)
env.SConscript('demos/SConscript', exports='env build_dir cache_dir')
env.SConscript('src/SConscript', exports='env build_dir cache_dir')

# Export files

if GetOption('production'):
    Default('html')

if 'html' in COMMAND_LINE_TARGETS or GetOption('production'):
    env.Execute('rm -rf $BASE_DIR/html')
env.Command('html', '$BUILD_DIR/index2.html',
            ['cp -r $BUILD_DIR $TARGET'])

def when_done():
    from SCons.Script import GetBuildFailures
    if not list(GetBuildFailures()):
        print("")
        print("Build successful!  Open or reload")
        print("     http://localhost:8081/")
        print("in your browser to see the result.")

atexit.register(when_done)


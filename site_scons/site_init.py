
import os

from subprocess import check_output


def cat_files(target, sources, delim=b''):
    with open(target, 'wb') as outfobj:
        for infile in sources:
            with open(str(infile), 'br') as infobj:
                outfobj.write(infobj.read())
                outfobj.write(delim)

def cat_css(target, source, env, for_signature=None):
    target = str(target[0])
    cat_files(target, source)

def cat_js(target, source, env, for_signature=None):
    target = str(target[0])
    cat_files(target, source, b';\n')

def minify(target, source, env, for_signature=None):
    target = str(target[0])
    src = str(source[0])
    if not env['MINIFY']:
        minimized = source[0].get_contents()
    elif src.endswith('.js'):
        minimized = check_output([env['UGLIFYJS'], '-m', '--', src])
    elif src.endswith('.css'):
        minimized = check_output([env['CLEANCSS'], '--skip-rebase', src])
    else:
        raise Exception("Don't know how to minify %s!" % src)
    if not os.path.exists(os.path.dirname(target)):
        os.makedirs(os.path.dirname(target))
    with open(target, 'wb') as fobj:
        fobj.write(minimized)


def TOOL_ADD_CAT(env):
    env['BUILDERS']['CatCSS'] = Builder(
        action=Action(cat_css, '$TARGET <-- $SOURCES'),
        suffix='.css', src_suffix='.css')
    env['BUILDERS']['CatJS'] = Builder(
        action=Action(cat_js, '$TARGET <-- $SOURCES'),
        suffix='.js', src_suffix='.js')
    env['BUILDERS']['Minify'] = Builder(
        action=Action(minify,
                      'minify $TARGET <-- $SOURCE'
                      if env['MINIFY'] else
                      'copy $TARGET <-- $SOURCE'))

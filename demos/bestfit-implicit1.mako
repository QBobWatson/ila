## -*- coffee -*-

<%inherit file="base2.mako"/>

<%block name="title">Implicit Function of Best Fit</%block>

<%block name="inline_style">
${parent.inline_style()}
  #eqn-here {
    color: var(--palette-red);
  }
</%block>

##

range = urlParams.get 'range', 'float', 10
rangeZ = urlParams.get 'rangez', 'float', range

# urlParams.func has to be of the form A*blah1(x)+B*blah2(x)-C*blah3(x)+D
# spaces are replaced by '+' (the reverse happens in urldecode)
funcStr = urlParams.func ? 'C*x+D'
funcStr = funcStr.replace /\s+/g, '+'
func = exprEval.Parser.parse funcStr
# x and y are function variables; the rest are parameters
params = []
zeroParams = {}
uniforms = {}
vars = ['x', 'y']
for letter in func.variables().sort()
    if letter in ['x', 'y']
        continue
    params.push letter
    zeroParams[letter] = 0
    uniforms[letter] =
        type:  'f'
        value: 0.0
numParams = params.length

# unit coordinate vectors in the parameters
units = []
for i in [0...numParams]
    obj = {}
    for param in params
        obj[param] = 0
    obj[params[i]] = 1
    units.push obj

# Target vectors
targets = []
i = 1
while urlParams["v#{i}"]?
    target = urlParams.get "v#{i}", 'float[]'
    target[2] = 0
    targets.push target
    i++
numTargets = targets.length

# Set up the linear equations
# The matrix has numParams columns, each with numTargets rows
matrix = ((0 for [0...numTargets]) for [0...numParams])
bvec = (0 for [0...numTargets])
xhat = (0 for [0...numParams])
bestfit = (x, y) -> 0
bestFitStr = ''
paramVals = {}

# Pass the function (almost) directly from the URL into the GPU
funcStrUrl = funcStr.replace /([a-zA-Z_]+)\^(\d+)/g, (match, p1, p2) ->
    (p1 for [0...parseInt(p2)]).join("*")
funcFragment = ''
for letter in params
    funcFragment += "uniform float #{letter};\n"
funcFragment += """\n
    float func(float x, float y) {
        return #{funcStrUrl};
    }
"""
# This should work for most TeX functions
for op of func.unaryOps
    if op.match /^[a-zA-Z]+$/
        funcStr = funcStr.replace(new RegExp(op, 'g'), "\\#{op}")

dot = (v1, v2) ->
    ret = 0
    for i in [0...v1.length]
        ret += v1[i] * v2[i]
    ret

updateCaption = () ->

solve = () ->
    # First have to figure out the coefficients of matrix and bvec
    for eqno in [0...numTargets]
        target = targets[eqno]
        linear = func.simplify
            x: target[0]
            y: target[1]
        # linear is now an (affine linear) function of the parameters only
        # First get the constant term
        constant = linear.evaluate zeroParams
        # Now get the linear terms
        for i in [0...numParams]
            matrix[i][eqno] = linear.evaluate(units[i]) - constant
        # The last coordinate of the target is the right-hand side of the equation
        bvec[eqno] = target[2] - constant
    # Now least-squares solve Ax=b
    ATA = ((dot(matrix[i], matrix[j]) for i in [0...numParams]) for j in [0...numParams])
    ATb = (dot(matrix[i], bvec) for i in [0...numParams])
    solver = rowReduce(ATA)[3]
    solver ATb, xhat
    # Substitute the parameters to get the best-fit function
    for letter, i in params
        paramVals[letter] = xhat[i]
        uniforms[letter].value = xhat[i]
    bestfit = func.simplify(paramVals).toJSFunction(vars.join ',')
    makeString()
    updateCaption()

makeString = () ->
    # Make a TeX string out of the function
    bestFitStr = funcStr
    for letter, i in params
        val = xhat[i]
        if val >= 0
            valAlone = val.toFixed 2
            valPlus  = "+#{valAlone}"
            valMinus = "-#{valAlone}"
        if val < 0
            val = -val
            valAlone = val.toFixed 2
            valPlus  = "-#{valAlone}"
            valMinus = "+#{valAlone}"
            valAlone = "-#{valAlone}"
        bestFitStr = bestFitStr.replace(
            new RegExp("\\+#{letter}\\*", 'g'), valPlus + '\\,')
        bestFitStr = bestFitStr.replace(
            new RegExp("\\+#{letter}", 'g'), valPlus)
        bestFitStr = bestFitStr.replace(
            new RegExp("-#{letter}\\*", 'g'), valMinus + '\\,')
        bestFitStr = bestFitStr.replace(
            new RegExp("-#{letter}", 'g'), valMinus)
        bestFitStr = bestFitStr.replace(
            new RegExp("#{letter}\\*", 'g'), valAlone + '\\,')
        bestFitStr = bestFitStr.replace(
            new RegExp("#{letter}", 'g'), valAlone)
        bestFitStr = bestFitStr.replace(/\*/g, '')

solve()


clipShader = \
    """
    // Enable STPQ mapping
    #define POSITION_STPQ
    void getPosition(inout vec4 xyzw, inout vec4 stpq) {
      // Store XYZ per vertex in STPQ
    stpq = xyzw;
    }
    """

# A point is drawn on the line if the function changes sign nearby.
# "nearby" means "on one of the points in toSample"
# glsl2 can't do statically initialized arrays as far as I can tell
numSamples = 8
radius = 0.02
samples = ''
for i in [0...numSamples]
    c = Math.cos(2*π*i/numSamples) * radius
    s = Math.sin(2*π*i/numSamples) * radius
    samples += "if(func(stpq.x + #{c.toFixed 8}, stpq.y + #{s.toFixed 8}) * val < 0.0)\n"
    samples += "    return rgba;\n"

curveFragment = \
    """
    // Enable STPQ mapping
    #define POSITION_STPQ

    vec4 getColor(vec4 rgba, inout vec4 stpq) {
        float val = func(stpq.x, stpq.y);

        #{samples}

        discard;
    }
    """


window.demo = new Demo2D {}, () ->
    window.mathbox = @mathbox

    view = @view
        axes:      true
        axisLabels: false
        grid:      true
        viewRange: [[-range,range],[-range,range]]

    ##################################################
    # (Unlabeled) points
    @labeledPoints view,
        name:      'targets'
        points:    targets
        colors:    (new Color("blue") for [0...numTargets])
        live:      true
        pointOpts: zIndex: 2

    ##################################################
    # Implicit curve
    # We don't compute the curve, so much as run a surface through a custom clip
    # shader, and only draw values near z=0.
    curve = view
        .shader code: clipShader
        .vertex pass: 'data'
        .shader
            code: funcFragment + "\n" + curveFragment
            uniforms: uniforms
        .fragment()
    curve
        .matrix
            channels: 2
            width:    2
            height:   2
            data:     [[[-range,-range], [-range,range]],
                       [[ range,-range], [ range,range]]]
        .surface
            color:   new Color("red").arr()
            opacity: 1.0
            fill:    true

    ##################################################
    # Dragging
    @draggable view,
        points:   targets
        postDrag: solve

    ##################################################
    # Caption
    str = '<p>Best-fit equation: <span id="eqn-here"></span></p>'
    @caption str

    bestFitElt = document.getElementById 'eqn-here'

    updateCaption = () =>
        katex.render "\\quad 0 = #{bestFitStr}", bestFitElt

    updateCaption()

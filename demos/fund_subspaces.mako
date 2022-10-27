## -*- coffee -*-

<%inherit file="base_diptych.mako"/>

<%block name="title">The Four Fundamental Subspaces</%block>

<%block name="inline_style">
${parent.inline_style()}
#matrix-here {
    text-align: center;
}
.overlay-text > p:last-child {
    text-align: center;
}
.mathbox-label p {
    text-align: center
}
</%block>

<%block name="label1">
<div class="mathbox-label">
<p>Row Picture</p>
<p id="row-picture-here"></p>
</div>
</%block>

<%block name="label2">
<div class="mathbox-label">
<p>Column Picture<p>
<p id="col-picture-here"></p>
</div>
</%block>

<%block name="overlay_text">
<div class="overlay-text">
  <p id="matrix-here"><span id="the-matrix"></span></p>
  <p>[Click and drag the rows and columns]</p>
</div>
</%block>


##################################################
# Globals

col1 = urlParams.get 'v1', 'float[]', [5, 3, -2]
col2 = urlParams.get 'v2', 'float[]', [3, -4, 1]
col3 = urlParams.get 'v3', 'float[]', [-1, 1, 7]

numRows = col1.length

if urlParams.v1?
    numCols = 1
if urlParams.v2?
    numCols = 2
if urlParams.v3?
    numCols = 3

rows = [[0,0,0][0...numCols],
        [0,0,0][0...numCols],
        [0,0,0][0...numCols]][0...numRows]
cols = [col1, col2, col3][0...numCols]

updateRows = () ->
    rows[0][0] = cols[0][0]
    if numCols > 1
        rows[0][1] = cols[1][0]
    if numCols > 2
        rows[0][2] = cols[2][0]
    if numRows > 1
        rows[1][0] = cols[0][1]
        if numCols > 1
            rows[1][1] = cols[1][1]
        if numCols > 2
            rows[1][2] = cols[2][1]
    if numRows > 2
        rows[2][0] = cols[0][2]
        if numCols > 1
            rows[2][1] = cols[1][2]
        if numCols > 2
            rows[2][2] = cols[2][2]
updateRows()

updateCols = () ->
    cols[0][0] = rows[0][0]
    if numRows > 1
        cols[0][1] = rows[1][0]
    if numRows > 2
        cols[0][2] = rows[2][0]
    if numCols > 1
        cols[1][0] = rows[0][1]
        if numRows > 1
            cols[1][1] = rows[1][1]
        if numRows > 2
            cols[1][2] = rows[2][1]
    if numCols > 2
        cols[2][0] = rows[0][2]
        if numRows > 1
            cols[2][1] = rows[1][2]
        if numRows > 2
            cols[2][2] = rows[2][2]


window.demo1 = new (if numCols == 3 then Demo else Demo2D) {
    mathbox: element: document.getElementById "mathbox1"
    scaleUI: true
}, () ->
    window.mathbox1 = @mathbox
    is2D = (numCols == 2)

    ##################################################
    # view
    @range = @urlParams.get 'range1', 'float', 10
    r = @range
    view = @view
        name:       'view1'
        viewRange:  [[-r,r], [-r,r], [-r,r]][0...numRows]
        grid:       false
        axes:       false

    @labels = ['r1', 'r2', 'r3'][0...numRows]
    @colors = [new Color("blue"), new Color("green"), new Color("brown")][0...numRows]

    ##################################################
    # Clip cube
    clipCube = @clipCube view,
        draw:   true
        hilite: not is2D
        material: new THREE.MeshBasicMaterial
            opacity:     0.5
            transparent: true
            visible:     false
            depthWrite:  false
            depthTest:   true

    unless is2D
        clipCube.installMesh()

    ##################################################
    # labeled vectors
    lVectors = rows.slice()
    lColors  = @colors.slice()
    lLabels  = @labels.slice()

    @labeledVectors view,
        vectors:       lVectors
        colors:        lColors
        labels:        lLabels
        live:          true
        zeroPoints:    true
        zeroThreshold: 0.05
        vectorOpts:    zIndex: 2
        labelOpts:     zIndex: 3
        zeroOpts:      zIndex: 3

    ##################################################
    # Subspace
    snapThreshold = 1.0 * @range / 10.0
    zeroThreshold = 0.00001

    @rowspaceColor = new Color "orange"
    @nullspaceColor = new Color "pink"

    subspace = @subspace
        vectors:       rows
        zeroThreshold: zeroThreshold
        live:          true
        range:         @range
        color:         @rowspaceColor
        # Lines before planes for transparency
        lineOpts:
            zOrder: 0
        surfaceOpts:
            zOrder: 1
    subspace.draw clipCube.clipped

    complement = @subspace
        name:          'complement'
        vectors:       subspace.complementFull(is2D)
        zeroThreshold: zeroThreshold
        live:          true
        range:         @range
        color:         @nullspaceColor
        pointOpts:     {size: 20, zIndex: 4}
        # Lines before planes for transparency
        lineOpts:
            zOrder: 0
        surfaceOpts:
            zOrder: 1
    complement.draw clipCube.clipped

    updateMesh = () =>
        return if is2D
        mesh = clipCube.mesh
        if complement?.dim == 3
            mesh.material.color = @nullspaceColor
            mesh.material.visible = true
        else if subspace.dim == 3
            mesh.material.color = @rowspaceColor.three()
            mesh.material.visible = true
        else
            mesh.material.visible = false
    updateMesh()

    ##################################################
    # Snap to subspace
    snapped = new THREE.Vector3()
    diff = new THREE.Vector3()
    ss0 = @subspace vectors: [[0, 0, 0]]
    ss1 = @subspace vectors: [[0, 0, 0]]
    ss2 = @subspace vectors: [[0, 0, 0], [0, 0, 0]]
    ss3 = @subspace vectors: [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    subspaces = [ss0, ss1, ss2, ss3]

    snap = (vec, vecs) =>
        ss = subspaces[vecs.length]
        if vecs.length > 0
            ss.setVecs vecs
        ss.project vec, snapped
        diff.copy(vec).sub snapped
        if diff.lengthSq() <= snapThreshold
            vec.copy snapped
            return true
        return false

    onDrag = (vec) ->
        updateCols()
        # Try snapping to 0
        return if snap vec, []
        indices = [0...numRows].filter (x) => x != @dragging
        others = (rows[i] for i in indices)
        # Try snapping to one of the other vectors
        for other in others
            return if snap vec, [other]
        # Try snapping to the span of the other vectors
        snap vec, others

    @postDrag = () =>
        subspace.setVecs rows
        complement?.setVecs subspace.complementFull(is2D)
        updateMesh()

    # Make the vectors draggable
    @draggable view,
        points: rows
        onDrag: onDrag
        postDrag: () =>
            @postDrag()
            window.demo2.postDrag()


window.demo2 = new (if numRows == 3 then Demo else Demo2D) {
    mathbox: element: document.getElementById "mathbox2"
    scaleUI: true
}, () ->
    window.mathbox2 = @mathbox
    is2D = (numRows == 2)

    ##################################################
    # view, axes
    @range = @urlParams.get 'range2', 'float', 10
    r = @range
    view = @view
        name:       'view2'
        viewRange:  [[-r,r], [-r,r], [-r,r]][0...numRows]
        grid:       false
        axes:       false

    @labels = ['c1', 'c2', 'c3'][0...numCols]
    @colors = [new Color("blue"), new Color("green"), new Color("brown")][0...numCols]

    ##################################################
    # Clip cube
    clipCube = @clipCube view,
        draw:   true
        hilite: not is2D
        material: new THREE.MeshBasicMaterial
            opacity:     0.5
            transparent: true
            visible:     false
            depthWrite:  false
            depthTest:   true

    unless is2D
        clipCube.installMesh()

    ##################################################
    # labeled vectors
    lVectors = cols.slice()
    lColors  = @colors.slice()
    lLabels  = @labels.slice()

    @labeledVectors view,
        vectors:       lVectors
        colors:        lColors
        labels:        lLabels
        live:          true
        zeroPoints:    true
        zeroThreshold: 0.05
        vectorOpts:    zIndex: 2
        labelOpts:     zIndex: 3
        zeroOpts:      zIndex: 3

    ##################################################
    # Subspace
    snapThreshold = 1.0 * @range / 10.0
    zeroThreshold = 0.00001

    @colspaceColor = new Color "violet"
    @leftnullspaceColor = new Color "red"

    subspace = @subspace
        vectors:       cols
        zeroThreshold: zeroThreshold
        live:          true
        range:         @range
        color:         @colspaceColor
        # Lines before planes for transparency
        lineOpts:
            zOrder: 0
        surfaceOpts:
            zOrder: 1
    subspace.draw clipCube.clipped

    complement = @subspace
        name:          'complement'
        vectors:       subspace.complementFull(is2D)
        zeroThreshold: zeroThreshold
        live:          true
        range:         @range
        color:         @leftnullspaceColor
        pointOpts:     {size: 20, zIndex: 4}
        # Lines before planes for transparency
        lineOpts:
            zOrder: 0
        surfaceOpts:
            zOrder: 1
    complement.draw clipCube.clipped

    updateMesh = () =>
        return if is2D
        mesh = clipCube.mesh
        if complement?.dim == 3
            mesh.material.color = @leftnullspaceColor
            mesh.material.visible = true
        else if subspace.dim == 3
            mesh.material.color = @colspaceColor.three()
            mesh.material.visible = true
        else
            mesh.material.visible = false
    updateMesh()

    ##################################################
    # Snap to subspace
    snapped = new THREE.Vector3()
    diff = new THREE.Vector3()
    ss0 = @subspace vectors: [[0, 0, 0]]
    ss1 = @subspace vectors: [[0, 0, 0]]
    ss2 = @subspace vectors: [[0, 0, 0], [0, 0, 0]]
    ss3 = @subspace vectors: [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    subspaces = [ss0, ss1, ss2, ss3]

    snap = (vec, vecs) =>
        ss = subspaces[vecs.length]
        if vecs.length > 0
            ss.setVecs vecs
        ss.project vec, snapped
        diff.copy(vec).sub snapped
        if diff.lengthSq() <= snapThreshold
            vec.copy snapped
            return true
        return false

    onDrag = (vec) ->
        updateRows()
        # Try snapping to 0
        return if snap vec, []
        indices = [0...numCols].filter (x) => x != @dragging
        others = (cols[i] for i in indices)
        # Try snapping to one of the other vectors
        for other in others
            return if snap vec, [other]
        # Try snapping to the span of the other vectors
        snap vec, others

    @postDrag = () =>
        subspace.setVecs cols
        complement?.setVecs subspace.complementFull(is2D)
        updateMesh()
        @updateCaption()

    # Make the vectors draggable
    @draggable view,
        points: cols
        onDrag: onDrag
        postDrag: () =>
            @postDrag()
            demo1.postDrag()

    @updateCaption = () ->
        str = @texMatrix cols,
            precision: 2
        katex.render 'A=' + str, document.getElementById('the-matrix')
    @updateCaption()

    katex.render """
          \\textcolor{#{demo1.rowspaceColor.str()}}{\\mathrm{Row}(A)} \\quad
          \\textcolor{#{demo1.nullspaceColor.str()}}{\\mathrm{Nul}(A)}
        """, document.getElementById('row-picture-here')
    katex.render """
          \\textcolor{#{@colspaceColor.str()}}{\\mathrm{Col}(A)} \\quad
          \\textcolor{#{@leftnullspaceColor.str()}}{\\mathrm{Nul}(A^T)}
        """, document.getElementById('col-picture-here')


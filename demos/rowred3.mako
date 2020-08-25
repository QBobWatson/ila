## -*- coffee -*-

<%inherit file="base_slideshow.mako"/>

<%block name="title">Row Reducing a Matrix</%block>
<%block name="first_caption">Let's row reduce this matrix.</%block>

## '

rrmat.setMatrix [[0, -7, -4,  2],
                 [2,  4,  6, 12],
                 [3,  1, -1, -2]]

window.slideshow =
    rrmat.slideshow()

    .caption "<b>Step 1a:</b> These are the nonzero entries in the first column."
    .setStyle(blink "red", [[1,0],[2,0]])
    .caption "<b>Step 1a:</b> These are the nonzero entries in the first column." \
             + " We need to move one of them to the first row."
    .break()

    .setStyle
        color:    "black"
        entries:  [[2,0]]
        duration: 0.2
    .rowSwap 0, 1
    .break()

    .caption "<b>Step 1b:</b> We kill this entry..."
    .setStyle(blink "blue", [[2,0]], 1)
    .caption "<b>Step 1b:</b>" \
             + " We kill this entry by subtracting 3/2 times the first row."
    .rowRep 0, -3/2, 2
    .setStyle
        color:    "black",
        entries:  [[2,0]]
        duration: 0.2
    .break()

    .caption "Let's get rid of this distracting line."
    .unAugment()
    .break()

    .caption """Now that the first column is clear except for the pivot in
                the first row, we can ignore the first row and the first
                column and concentrate on the rest of the matrix."""
    .setStyle
        color:    "rgb(200,200,200)",
        entries:  [[0,0],[0,1],[0,2],[0,3],[1,0],[2,0]]
        duration: 0.5
    .break()

    .caption "<b>Step 2b:</b> To kill this entry..."
    .setStyle(blink "blue", [[2,1]], 1)
    .caption """<b>Step 2b:</b>
                To kill this entry, we subtract 5/7 times row 2 from row 3."""
    .rowRep 1, -5/7, 2
    .setStyle
        color:    "black"
        entries:  [[2,1]]
        duration: 0.2
    .break()

    .caption """Now the second column (the part we care about) is clear
                except for the pivot, so we would ignore the
                second column and second row and continue recursively, if
                there were more rows and columns."""
    .setStyle
        color:    "rgb(200,200,200)"
        entries:  [[1,1],[1,2],[1,3],[2,1]]
        duration: 1
    .break()

    .caption "Notice that the matrix is now in <i>row echelon form</i>."
    .setStyle [{
        color:    "red"
        entries:  [[0,0],[1,1],[2,2]]
        duration: 1
        }, {
        color:    "black"
        entries:  [[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]]
        duration: 1
        }]
    .break()

    .caption """To put the matrix in <i>reduced</i> row echelon form, we
                need to make the <span style=\"color: red;\">pivots</span>
               equal to one"""
    .setStyle(blink "red", [[0,0],[1,1],[2,2]], 1)
    .break()

    .caption """To put the matrix in <i>reduced</i> row echelon form, we
                need to make the <span style=\"color: red;\">pivots</span>
                equal to one, then kill the <span style=\"color: blue;\">entries</span>
                above the pivots."""
    .setStyle(blink "blue", [[0,1],[0,2],[1,2]], 1)
    .break()

    .caption """We start with the last pivot"""
    .setStyle [
        color:   "black"
        entries: [[0,0],[0,1],[0,2],[1,1],[1,2]]
        ].concat(blink "red", [[2, 2]], 1)

    .caption """We start with the last pivot, and set it equal to 1 by scaling."""
    .rowMult 2, -7/50
    .break()

    .caption "Now we clear the third column..."
    .setStyle(blink "blue", [[0,2],[1,2]], 1)
    .caption """Now we clear the third column column using row replacement."""
    .rowRep 2, 4, 1
    .rowRep 2, -6, 0
    .break()

    .caption """We scale the next pivot to make it equal to 1"""
    .setStyle [
        color:   "black"
        entries: [[0,2],[1,2],[2,2]]
        ].concat(blink "red", [[1,1]], 1)
    .rowMult 1, -1/7
    .break()

    .caption "Then we clear the second column..."
    .setStyle(blink "blue", [[0,1]], 1)
    .caption "Then we clear the second column using row replacement."
    .rowRep 1, -4, 0
    .break()

    .caption """Finally we scale the first pivot to make it equal to 1"""
    .setStyle [
        color:   "black"
        entries: [[0,1],[1,1]]
        ].concat(blink "red", [[0,0]], 1)
    .rowMult 0, 1/2
    .break()

    .caption "The matrix is now in reduced row echelon form!"
    .setStyle [{
        color:    "black"
        entries:  [[1,0],[2,0],[2,1],[0,1]]
        duration: 1}, {
        color: "red"
        entries: [[0,0],[1,1],[2,2]]
        duration: 1
        }]
    .break()

    .caption "Add back the divider..."
    .reAugment()
    .caption "Add back the divider, and we're done!"
    .setStyle
        transform: "rotate(360deg)"
        entries:   [[0,0],[0,1],[0,2],[0,3],
                    [1,0],[1,1],[1,2],[1,3],
                    [2,0],[2,1],[2,2],[2,3]]
        duration:  1.5
        timing:    'linear'
    .break()

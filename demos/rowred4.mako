
<%inherit file="base_slideshow.mako"/>

<%block name="title">Gaussian Elimination</%block>
<%block name="first_caption">Let's perform Gaussian elimination on this matrix.</%block>

## '

rrmat = window.rrmat = new RRMatrix 4, 4, view, mathbox,
    augmentCol:     2
    startAugmented: true

rrmat.setMatrix [
    [ 1,   7,  2,   4],
    [-2, -14, -4, -14],
    [ 2,  18,  2,   8],
    [-1,  -5, -3,  -7]
]

window.slideshow =
    rrmat.slideshow()

    .caption "<b>Step 1b:</b> We kill this entry..."
    .setStyle(blink "blue", [[1,0]], 1)
    .caption "<b>Step 1b:</b>" \
             + " We kill this entry by adding 2 times the first row."
    .rowRep 0, 2, 1
    .setStyle
        color:    "black",
        entries:  [[1,0]]
        duration: 0.2
    .break()

    .caption "<b>Step 1b:</b> We clear the rest of the first column..."
    .setStyle(blink "blue", [[2,0],[3,0]], 1)
    .caption "<b>Step 1b:</b>" \
             + " We clear the rest of the first column in the same way."
    .rowRep 0, -2, 2
    .rowRep 0, 1, 3
    .setStyle
        color:    "black",
        entries:  [[2,0],[3,0]]
        duration: 0.2
    .break()

    .caption """Now that the first column is clear except for the pivot in
                the first row, we can ignore the first row and the first
                column and concentrate on the rest of the matrix."""
    .setStyle
        color:    "rgb(200,200,200)",
        entries:  [[0,0],[0,1],[0,2],[0,3],[1,0],[2,0],[3,0]]
        duration: 0.5
    .break()

    .caption "Let's get rid of this distracting line."
    .unAugment()
    .break()

    .caption "<b>Step 2a:</b> These are the nonzero entries in the second column."
    .setStyle(blink "red", [[2,1],[3,1]])
    .caption "<b>Step 2a:</b> These are the nonzero entries in the second column." \
             + " We need to move one of them to the second row."
    .break()

    .setStyle
        color:    "black"
        entries:  [[2,1],[3,1]]
        duration: 0.2
    .rowSwap 1, 3
    .break()

    .caption "<b>Step 2b:</b> We kill this entry..."
    .setStyle(blink "blue", [[2,1]], 1)
    .caption "<b>Step 1b:</b>" \
             + " We kill this entry by subtracting 2 times the second row."
    .rowRep 1, -2, 2
    .setStyle
        color:    "black",
        entries:  [[2,1]]
        duration: 0.2
    .break()

    .caption "Now we have completed steps 1 and 2."
    .setStyle
        color:    "rgb(200,200,200)",
        entries:  [[1,1],[1,2],[1,3],[2,1],[2,2],[3,1],[3,2]]
        duration: 0.5
    .break()
    .caption """Now we have completed steps 1 and 2.  There is no pivot in the third
                column; our third pivot is in the fourth column."""
    .setStyle(blink "red", [[2,3]])
    .break()

    .setStyle
        color:    "black",
        entries:  [[2,3]]
        duration: 0.2
    .caption "<b>Step 2b:</b> We kill this entry..."
    .setStyle(blink "blue", [[3,3]], 1)
    .caption "<b>Step 1b:</b>" \
             + " We kill this entry by adding the third row."
    .rowRep 2, 1, 3
    .setStyle
        color:    "black",
        entries:  [[3,3]]
        duration: 0.2
    .break()

    .caption "Now our matrix is in REF!"
    .setStyle
        color:    "black"
        entries:  [[0,0],[0,1],[0,2],[0,3],[1,0],[2,0],[3,0],
                   [1,1],[1,2],[1,3],[2,1],[2,2],[3,1],[3,2]]
        duration: 0.2
    .setStyle
        color:    "red"
        entries:  [[0,0],[1,1],[2,3]]
        duration: 0.2
    .reAugment()
    .break()
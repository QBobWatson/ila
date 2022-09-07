
<%inherit file="base_slideshow.mako"/>

<%block name="title">Jordan Substitution</%block>
<%block name="first_caption">Let's perform Jordan substitution on this matrix.</%block>

## '

rrmat = window.rrmat = new RRMatrix 3, 4, view, mathbox,
    augmentCol:     2
    startAugmented: true

rrmat.setMatrix [
    [1,  2,   3,   6],
    [0, -5, -10, -20],
    [0,  0,  10,  30]
]

window.slideshow =
    rrmat.slideshow()

    .caption "<b>Step 1:</b> To make this pivot into a 1..."
    .setStyle(blink "red", [[2,2]], 1)
    .caption "<b>Step 1:</b>" \
             + " To make this pivot into a 1, we divide the third row by 10."
    .rowMult 2, 1/10
    .break()

    .caption "<b>Step 2:</b> We clear these entries..."
    .setStyle(blink "blue", [[1,2],[0,2]], 1)
    .caption "<b>Step 2:</b> We clear these entries by adding multiples of the third row."
    .rowRep 2, -3, 0
    .rowRep 2, 10, 1
    .setStyle
        color:    "black",
        entries:  [[0,2],[1,2],[2,2]]
        duration: 0.2
    .break()

    .caption "<b>Step 1:</b> To make this pivot into a 1..."
    .setStyle(blink "red", [[1,1]], 1)
    .caption "<b>Step 1:</b>" \
             + " To make this pivot into a 1, we divide the second row by -5."
    .rowMult 1, -1/5
    .break()

    .caption "<b>Step 2:</b> We clear this entry..."
    .setStyle(blink "blue", [[0,1]], 1)
    .caption "<b>Step 2:</b> We clear this entry by subtracting 2 times the second row."
    .rowRep 1, -2, 0
    .setStyle
        color:    "black",
        entries:  [[0,1],[1,1]]
        duration: 0.2
    .break()

    .caption "This matrix is in reduced row echelon form!"
    .setStyle
        transform: "rotate(360deg)"
        entries: [[0,0],[0,1],[0,2],[0,3],
                  [1,0],[1,1],[1,2],[1,3],
                  [2,0],[2,1],[2,2],[2,3]]
        duration: 1.5
        timing: 'linear'
    .break()

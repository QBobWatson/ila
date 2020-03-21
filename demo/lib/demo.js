'use strict';

import OrbitControls from './orbitcontrols.js';

////////////////////////////////////////////////////////////////////////////////
// * Color

// Keep this in sync with jdr-tikz.sty
let palette = {
    red:    [0.8941, 0.1020, 0.1098],
    blue:   [0.2157, 0.4941, 0.7216],
    green:  [0.3020, 0.6863, 0.2902],
    violet: [0.5961, 0.3059, 0.6392],
    orange: [1.0000, 0.4980, 0.0000],
    yellow: [0.7000, 0.7000, 0.0000],
    brown:  [0.6510, 0.3373, 0.1569],
    pink:   [0.9686, 0.5059, 0.7490]
};

class Color {
    constructor(...args) {
        this.r = this.g = this.b = 1.0;
        if(args.length > 0)
            this.set(...args);
    }

    toArray(...args) {
        return [this.r, this.g, this.b, ...args];
    }

    get rgb() {
        return this.toArray();
    }

    set rgb(arr) {
        [this.r, this.g, this.b] = arr;
        this.r = Math.min(1.0, Math.max(0.0, this.r));
        this.g = Math.min(1.0, Math.max(0.0, this.g));
        this.b = Math.min(1.0, Math.max(0.0, this.b));
    }

    get three() {
        return new THREE.Color(this.r, this.g, this.b);
    }

    set three(color) {
        this.rgb = [color.r, color.g, color.b];
    }

    get num() {
        return (Math.round(this.r*255) << 16)
            |  (Math.round(this.g*255) << 8)
            |  (Math.round(this.b*255) << 0);
    }

    set num(number) {
        number = Math.floor(number);
        this.rgb = [
            (number >> 16 & 255) / 255,
            (number >> 8  & 255) / 255,
            (number       & 255) / 255
        ];
    }

    get str() {
        return '#' + ('000000' + this.hex.toString(16)).slice(-6);
    }

    set str(str) {
        let color;
        if(color = /^rgb\((\d+), ?(\d+), ?(\d+)\)$/i.exec(str))
            // rgb(255,0,0)
            this.rgb = [
                Math.min(255, parseInt(color[1], 10)) / 255,
                Math.min(255, parseInt(color[2], 10)) / 255,
                Math.min(255, parseInt(color[3], 10)) / 255
            ];
        else if(color = /^rgb\((\d+)\%, ?(\d+)\%, ?(\d+)\%\)$/i.exec(str))
            // rgb(100%, 0%, 0%)
            this.rgb = [
                Math.min(100, parseInt(color[1], 10)) / 100,
                Math.min(100, parseInt(color[2], 10)) / 100,
                Math.min(100, parseInt(color[3], 10)) / 100
            ];
        else if(color = /^\#([0-9a-f]{6})$/i.exec(str))
            // #ff0000
            this.num = parseInt(color[1], 16);
        else if(color = /^\#([0-9a-f])([0-9a-f])([0-9a-f])$/i.exec(str))
            // #f00
            this.num = parseInt(color[1]+color[1]+
                                color[2]+color[2]+
                                color[3]+color[3], 16);
        else if(/^(\w+)$/i.test(str))
            // red
            this.rgb = palette[str];
        else
            throw 'Invalid argument to Color.str';
    }

    get hsl() {
        let {r, g, b} = this;
        let max = Math.max(r, g, b);
        let min = Math.min(r, g, b);
        let l = (max + min) / 2;
        let h = 0.0, s = 0.0;
        if(max != min) {
            let d = max - min;
            s = l > 0.5 ? d / (2 - max - min) : d / (max + min);
            switch(max) {
            case r:
                h = (g - b) / d + (g < b ? 6 : 0);
                break;
            case g:
                h = (b - r) / d + 2;
                break;
            case b:
                h = (r - g) / d + 4;
                break;
            }
            h /= 6;
        }
        return [h, s, l];
    }

    set hsl(hsl) {
        let [h, s, l] = hsl;

        if(s == 0.0) {
            this.rgb = [1.0, 1.0, 1.0];
            return;
        }

        let hue2rgb = (p, q, t) => {
            if(t < 0) t += 1;
            if(t > 1) t -= 1;
            if(t < 1/6) return p + (q - p) * 6 * t;
            if(t < 1/2) return q;
            if(t < 2/3) return p + (q - p) * (2/3 - t) * 6;
            return p;
        };

        let q = l < 0.5 ? l * (1 + s) : l + s - l * s;
        let p = 2 * l - q;

        this.rgb = [
            hue2rgb(p, q, h + 1/3),
            hue2rgb(p, q, h),
            hue2rgb(p, q, h - 1/3)
        ];
    }

    set(...args) {
        // Copied from THREE.js
        let r = 1.0, g = 1.0, b = 1.0;
        if(args.length >= 3)
            this.rgb = args;
        else if(Array.isArray(args[0]))
            this.rgb = args[0];
        else if(args[0] instanceof Color || args[0] instanceof THREE.Color)
            this.color = args[0];
        else if(typeof(args[0]) === 'number')
            this.num = args[0];
        else if(typeof(args[0]) == 'string')
            this.str = args[0];
        else
            throw 'Invalid argument to Color constructor';
        return this;
    }

    toString() {
        return `Color(${this.r}, ${this.g}, ${this.b})`;
    }

    brighten(pct) {
        let [h, s, l] = this.hsl;
        let x = new Color();
        x.hsl = [h, s, l+pct];
        return x;
    }

    darken(pct) {
        return this.brighten(-pct);
    }
}


////////////////////////////////////////////////////////////////////////////////
// * View

// Wrapper for a mathbox cartesian view (2D or 3D).
// Options:
//     name: ids and classes are prefixed with `${name}-`
//     viewRange: range option for the view.  Determines number of dimensions.
//     viewScale: scale option for the view
//     doAxes: construct the axes
//     axisOpts: options to mathbox.axis
//     doGrid: construct a grid
//     gridOpts: options to mathbox.grid
//     axisLabels: draw axis labels (x, y, z)
//     labelOpts: options to mathbox.label
//     viewOpts: options to mathbox.cartesian

class View {
    constructor(mathbox, {
        name       = "view",
        viewRange  = [[-10, 10], [-10, 10], [-10, 10]],
        viewScale  = [1, 1, 1],
        doAxes     = true,
        axisOpts   = {},
        doGrid     = true,
        gridOpts   = {},
        axisLabels = true,
        labelOpts  = {},
        viewOpts   = {}
    }={}) {
        this.numDims = viewRange.length;

        this.view = mathbox.cartesian({
            range: viewRange,
            scale: viewScale,
            id:    `${name}-view`,
            ...viewOpts
        });

        if(doAxes) {
            axisOpts = {
                classes: [`${name}-axes`],
                end:     true,
                width:   3,
                depth:   1,
                color:   "black",
                opacity: 0.5,
                zBias:   -1,
                size:    5,
                ...axisOpts
            };
            for(let i = 1; i <= this.numDims; ++i) {
                axisOpts.axis = i;
                this.view.axis(axisOpts);
            }

            if(axisLabels)
                this.view
                    .array({
                        channels: this.numDims,
                        width:    this.numDims,
                        live:     false,
                        expr: (emit, i) => {
                            let arr = [];
                            for(let j = 0; j < this.numDims; ++j) {
                                if(i == j) arr.push(viewRange[i][1] * 1.04);
                                else       arr.push(0.0);
                            }
                            emit(...arr);
                        }
                    }).text({
                        live:  false,
                        width: this.numDims,
                        data:  ['x', 'y', 'z'].splice(0, this.numDims)
                    }).label({
                        classes:    [`${name}-axes`],
                        size:       20,
                        color:      "black",
                        opacity:    0.5,
                        outline:    0,
                        background: [0,0,0,0],
                        offset:     [0, 0],
                        ...labelOpts
                    });
        }

        if(doGrid)
            this.view.grid({
                classes: [`${name}-axes`, `${name}-grid`],
                axes:    [1, 2],
                width:   2,
                depth:   1,
                color:   "black",
                opacity: 0.25,
                zBias:   0,
                ...gridOpts
            });
    }
}


////////////////////////////////////////////////////////////////////////////////
// * Demo

// Class for mathbox-based demos.  Parameters:
//    mathboxOpts: passed to the mathBox constructor
//    cameraOpts: passed to mathbox.camera()
//    clearColor: THREE's clear color
//    clearOpacity: THREE's clear opacity
//    focusDist: mathbox focus distance
//    scaleUI: whether to scale focusDist by min(width, height)/1000
//    dims: 2 or 3
//    preload: images to preload

class Demo {
    constructor(callback, {
        mathboxOpts  = {},
        cameraOpts   = {},
        clearColor   = 0xffffff,
        clearOpacity = 1.0,
        focusDist    = 1.5,
        scaleUI      = true,
        dims         = 3,
        preload      = {}
    }={}) {

        mathboxOpts = {
            plugins: ['core', 'controls', 'cursor'],
            controls: {
                klass: OrbitControls,
                parameters: {
                    noKeys: true
                }
            },
            mathbox: {
                inspect: false
            },
            splash: {
                fancy: true,
                color: "blue"
            },
            ...mathboxOpts
        };

        cameraOpts = {
            proxy:    true,
            position: [3, 1.5, 1.5],
            lookAt:   [0, 0, 0],
            up:       [0, 0, 1],
            ...cameraOpts
        };

        let onPreloaded = () => {
            this.mathbox = mathBox(mathboxOpts);
            this.three = this.mathbox.three;
            this.three.renderer.setClearColor(
                new Color(clearColor).three, clearOpacity);
            this.camera = this.mathbox.camera(cameraOpts)[0].controller.camera;
            this.controls = this.three.controls;
            this.controls && this.controls.updateCamera
                && this.controls.updateCamera();
            this.canvas = this.mathbox._context.canvas;
            if(scaleUI)
                this.mathbox.bind(
                    'focus', () => focusDist / 1000 *
                        Math.min(this.canvas.clientWidth, this.canvas.clientHeight));
            else
                this.mathbox.set('focus', focusDist);
            callback(this);
        };

        let toPreload = 0;
        for(let [key, val] of Object.entries(preload)) {
            toPreload++;
            let image = new Image();
            this[key] = image;
            image.src = val;
            image.addEventListener('load', () => {
                if(--toPreload == 0)
                    onPreloaded();
            });
        }
        toPreload || onPreloaded();
    }

    view(opts={}) {
        return new View(this.mathbox, opts);
    }

}


export { Demo, Color };

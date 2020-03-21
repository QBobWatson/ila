
// Modified from three.js OrbitControls:
//  * Can set 'up'
//  * Can control multiple cameras

export default class OrbitControls {
    constructor(camera, domElement=document) {
        THREE.EventDispatcher.prototype.apply(this);

        this.camera = camera;
        this.domElement = domElement;

        this.enabled         = true;
        this.target          = new THREE.Vector3();
        this.noZoom          = false;
        this.zoomSpeed       = 1.0;
        this.minDistance     = 0;
        this.maxDistance     = Infinity;
        this.noRotate        = false;
        this.rotateSpeed     = 1.0;
        this.noPan           = false;
        this.keyPanSpeed     = 7.0;
        this.autoRotate      = false;
        this.autoRotateSpeed = 2.0;
        this.minPolarAngle   = 0;
        this.maxPolarAngle   = Math.PI;
        this.noKeys          = true;

        this.keys = {
            LEFT:   37,
            UP:     38,
            RIGHT:  39,
            BOTTOM: 40
        };

        this.clones = [];

        // Internal state
        this.EPS = 0.000001;

        this.rotateStart = new THREE.Vector2();
        this.rotateEnd   = new THREE.Vector2();
        this.rotateDelta = new THREE.Vector2();

        this.panStart   = new THREE.Vector2();
        this.panEnd     = new THREE.Vector2();
        this.panDelta   = new THREE.Vector2();
        this.panOffset  = new THREE.Vector3();
        this.panCurrent = new THREE.Vector3();

        this.offset = new THREE.Vector3();

        this.dollyStart = new THREE.Vector2();
        this.dollyEnd   = new THREE.Vector2();
        this.dollyDelta = new THREE.Vector2();

        this.phiDelta   = 0;
        this.thetaDelta = 0;
        this.scale      = 1;

        this.lastPosition = new THREE.Vector3();

        this.STATE = {
            NONE:        -1,
            ROTATE:       0,
            DOLLY:        1,
            PAN:          2,
            TOUCH_ROTATE: 3,
            TOUCH_DOLLY:  4,
            TOUCH_PAN:    5
        };

        // for reset
        this.target0   = this.target.clone();
        this.position0 = this.camera.position.clone();

        this.updateCamera();

        // events
        this.changeEvent = {type: 'change'};
        this.startEvent  = {type: 'start'};
        this.endEvent    = {type: 'end'};

        // install listeners
        let de = this.domElement;
        de.addEventListener('contextmenu', event => event.preventDefault(), false);
        de.addEventListener('mousedown',   this.onMouseDown.bind(this),     false);
        de.addEventListener('mousewheel',  this.onMouseWheel.bind(this),    false);
        de.addEventListener('touchstart',  this.touchStart.bind(this),      false);
        window.addEventListener('keydown', this.onKeyDown.bind(this),       false);

        // force an update at start
        this.update();
    }

    enable(val) {
        this.enabled = val;
        if(this.enabled) return;
        let de = document.documentElement, s = this.STATE;
        if(this.state in [s.ROTATE, s.DOLLY, s.PAN]) {
            de.removeEventListener('mousemove', this.onMouseMove, false);
            de.removeEventListener('mouseup',   this.onMouseUp,   false);
            this.dispatchEvent(this.endEvent);
        } else if(this.state in [s.TOUCH_ROTATE, s.TOUCH_DOLLY, s.TOUCH_PAN]) {
            de.removeEventListener('touchend',    this.touchEnd,  false);
            de.removeEventListener('touchmove',   this.touchMove, false);
            de.removeEventListener('touchcancel', this.touchEnd,  false);
            this.dispatchEvent(this.endEvent);
        }
        this.state = s.NONE;
    }

    updateCamera() {
        // so camera.up is the orbit axis
        this.quat = new THREE.Quaternion().setFromUnitVectors(
            this.camera.up, new THREE.Vector3(0, 1, 0));
        this.quatInverse = this.quat.clone().inverse();
        this.update();
    }

    getAutoRotationAngle() {
        return 2 * Math.PI / 60 / 60 * this.autoRotateSpeed;
    }

    rotateLeft(angle=this.getAutoRotationAngle()) {
        this.thetaDelta -= angle;
    }

    rotateUp(angle=this.getAutoRotationAngle()) {
        this.phiDelta -= angle;
    }

    getZoomScale() {
        return Math.pow(0.95, this.zoomSpeed);
    }

    dollyIn(dollyScale=this.getZoomScale()) {
        this.scale /= dollyScale;
    }

    dollyOut(dollyScale=this.getZoomScale()) {
        this.scale *= dollyScale;
    }

    // pass in distance in world space to move left
    panLeft(distance) {
        let te = this.camera.matrix.elements
        // get X column of matrix
        this.panOffset.set(te[0], te[1], te[2]);
        this.panOffset.multiplyScalar(-distance);
        this.panCurrent.add(this.panOffset);
    }

    // pass in distance in world space to move up
    panUp(distance) {
        let te = this.camera.matrix.elements
        // get Y column of matrix
        this.panOffset.set(te[4], te[5], te[6]);
        this.panOffset.multiplyScalar(distance);
        this.panCurrent.add(this.panOffset);
    }

    // pass in x,y of change desired in pixel space,
    // right and down are positive
    pan(deltaX, deltaY) {
        let element = this.domElement == document ? document.body : this.domElement;
        let cam = this.camera;
        if(cam.fov !== undefined) {
            // perspective
            let position       = cam.position;
            let offset         = position.clone().sub(this.target);
            let targetDistance = offset.length();
            // half of the fov is center to top of screen
            targetDistance *= Math.tan(cam.fov/2 * Math.PI / 180.0);
            // we actually don't use screenWidth, since perspective camera
            // is fixed to screen height
            this.panLeft(2 * deltaX * targetDistance / element.clientHeight);
            this.panUp(2 * deltaY * targetDistance / element.clientHeight);
        } else if(cam.top !== undefined) {
            // orthographic
            this.panLeft(deltaX * (cam.right - cam.left)   / element.clientWidth)
            this.panUp  (deltaY * (cam.top -   cam.bottom) / element.clientHeight)
        } else
            console.warn('WARNING: OrbitControls encountered unknown camera type; pan disabled');
    }

    update(delta, state) {
        if(state == undefined)
            state = this;
        else {
            for(let clone of this.clones)
                clone.update(0, this);
        }

        let {thetaDelta, phiDelta, panCurrent, scale} = state;

        let position = this.camera.position;
        this.offset.copy(position).sub(this.target);
        // rotate offset to "y-axis-is-up" space
        this.offset.applyQuaternion(this.quat);
        let {x, y, z} = this.offset;
        // angle from z-axis around y-axis
        let theta = Math.atan2(x, z);
        // angle from y-axis
        let phi = Math.atan2(Math.sqrt(x * x + z * z), y);

        if(this.autoRotate) this.rotateLeft();

        theta += thetaDelta;
        phi += phiDelta;
        // restrict phi to be between desired limits
        phi = Math.max(this.minPolarAngle, Math.min(this.maxPolarAngle, phi));
        // restrict phi to be between EPS and PI-EPS
        phi = Math.max(this.EPS, Math.min(Math.PI - this.EPS, phi));
        let radius = this.offset.length() * scale;
        // restrict radius to be between desired limits
        radius = Math.max(this.minDistance, Math.min(this.maxDistance, radius));
        // move target to panned location
        this.target.add(panCurrent);
        this.offset.x = radius * Math.sin(phi) * Math.sin(theta);
        this.offset.y = radius * Math.cos(phi);
        this.offset.z = radius * Math.sin(phi) * Math.cos(theta);
        // rotate offset back to "camera-up-vector-is-up" space
        this.offset.applyQuaternion(this.quatInverse);
        // Update camera
        position.copy(this.target).add(this.offset);
        this.camera.lookAt(this.target);

        this.thetaDelta = 0;
        this.phiDelta = 0;
        this.scale = 1;
        this.panCurrent.set(0, 0, 0);

        if(this.lastPosition.distanceToSquared(position) > this.EPS)
            this.dispatchEvent(this.changeEvent);
        this.lastPosition.copy(position);
    }

    reset() {
        this.state = this.STATE.NONE;
        this.target.copy(this.target0);
        this.camera.position.copy(this.position0);
        this.update();
    }

    onMouseDown(event) {
        if(!this.enabled) return;
        event.preventDefault();

        switch(event.button) {
        case 0:
            if(this.noRotate) return;
            this.state = this.STATE.ROTATE;
            this.rotateStart.set(event.clientX, event.clientY);
            break;
        case 1:
            if(this.noZoom) return;
            this.state = this.STATE.DOLLY;
            this.dollyStart.set(event.clientX, event.clientY);
            break;
        case 2:
            if(this.noPan) return;
            this.state = this.STATE.PAN;
            this.panStart.set(event.clientX, event.clientY);
            break;
        default:
            return;
        }

        let de = document.documentElement;
        de.addEventListener('mousemove', this.onMouseMove.bind(this), false);
        de.addEventListener('mouseup',   this.onMouseUp.bind(this),   false);
        this.dispatchEvent(this.startEvent);
    }

    onMouseMove(event) {
        if(!this.enabled) return;
        event.preventDefault()

        let element = this.domElement == document ? document.body : this.domElement;

        switch(this.state) {
        case this.STATE.ROTATE:
            if(this.noRotate) return;

            this.rotateEnd.set(event.clientX, event.clientY);
            this.rotateDelta.subVectors(this.rotateEnd, this.rotateStart);
            // rotating across whole screen goes 360 degrees around
            this.rotateLeft(2 * Math.PI * this.rotateDelta.x
                            / element.clientWidth * this.rotateSpeed);
            // rotating up and down along whole screen attempts to go 360,
            // but limited to 180
            this.rotateUp(2 * Math.PI * this.rotateDelta.y
                          / element.clientHeight * this.rotateSpeed);
            this.rotateStart.copy(this.rotateEnd);
            break;

        case this.STATE.DOLLY:
            if(this.noZoom) return;

            this.dollyEnd.set(event.clientX, event.clientY);
            this.dollyDelta.subVectors(this.dollyEnd, this.dollyStart);
            if(this.dollyDelta.y > 0) this.dollyIn();
            else this.dollyOut();
            this.dollyStart.copy(this.dollyEnd);
            break;

        case this.STATE.PAN:
            if(this.noPan) return;

            this.panEnd.set(event.clientX, event.clientY);
            this.panDelta.subVectors(this.panEnd, this.panStart);
            this.pan(this.panDelta.x, this.panDelta.y);
            this.panStart.copy(this.panEnd);
            break;

        default:
            return;
        }

        this.update();
    }

    onMouseUp() {
        if(!this.enabled) return;

        let de = document.documentElement;
        de.removeEventListener('mousemove', this.onMouseMove, false);
        de.removeEventListener('mouseup',   this.onMouseUp,   false);
        this.dispatchEvent(this.endEvent);
        this.state = this.STATE.NONE;
    }

    onMouseWheel(event) {
        if(!this.enabled || this.noZoom) return;
        event.preventDefault();
        event.stopPropagation();

        let delta = event.wheelDelta !== undefined ? event.wheelDelta : -event.detail;
        if(delta > 0) this.dollyOut();
        else this.dollyIn();
        this.update();
        this.dispatchEvent(this.startEvent);
        this.dispatchEvent(this.endEvent);
    }

    onKeyDown(event) {
        if(!this.enabled || this.noKeys || this.noPan) return;

        switch(event.keyCode) {
        case this.keys.UP:      this.pan(0,  this.keyPanSpeed); break;
        case this.keys.BOTTOM:  this.pan(0, -this.keyPanSpeed); break;
        case this.keys.LEFT:    this.pan( this.keyPanSpeed, 0); break;
        case this.keys.RIGHT:   this.pan(-this.keyPanSpeed, 0); break;
        default: return;
        }

        this.update();
    }

    touchStart(event) {
        if(!this.enabled) return;
        event.preventDefault();

        switch(event.touches.length) {
        case 1:
            // one-fingered touch: rotate
            if(this.noRotate) return;
            this.state = this.STATE.TOUCH_ROTATE;
            this.rotateStart.set(event.touches[0].clientX, event.touches[0].clientY);
            break;
        case 2:
            // two-fingered touch: dolly
            if(this.noZoom) return;
            this.state = this.STATE.TOUCH_DOLLY;
            let dx = event.touches[0].clientX - event.touches[1].clientX;
            let dy = event.touches[0].clientY - event.touches[1].clientY;
            let distance = Math.sqrt(dx * dx + dy * dy);
            this.dollyStart.set(0, distance);
            break;
        case 3:
            // three-fingered touch: pan
            if(this.noPan) return;
            this.state = this.STATE.TOUCH_PAN;
            this.panStart.set(event.touches[0].clientX, event.touches[0].clientY);
            break;
        default:
            return;
        }

        let de = document.documentElement;
        de.addEventListener('touchend',    this.touchEnd.bind(this),  false);
        de.addEventListener('touchmove',   this.touchMove.bind(this), false);
        de.addEventListener('touchcancel', this.touchEnd.bind(this),  false);
        this.dispatchEvent(this.startEvent);
    }

    touchMove(event) {
        if(!this.enabled) return;
        event.preventDefault();
        event.stopPropagation();

        let element = this.domElement == document ? document.body : this.domElement;

        switch(event.touches.length) {
        case 1:
            // one-fingered touch: rotate
            if(this.noRotate || this.state != this.STATE.TOUCH_ROTATE) return;
            this.rotateEnd.set(event.touches[0].clientX, event.touches[0].clientY);
            this.rotateDelta.subVectors(this.rotateEnd, this.rotateStart);
            // rotating across whole screen goes 360 degrees around
            this.rotateLeft(2 * Math.PI * this.rotateDelta.x
                            / element.clientWidth * this.rotateSpeed);
            // rotating up and down along whole screen attempts to go 360,
            // but limited to 180
            this.rotateUp(2 * Math.PI * this.rotateDelta.y
                          / element.clientHeight * this.rotateSpeed);
            this.rotateStart.copy(this.rotateEnd);
            break;
        case 2:
            // two-fingered touch: dolly
            if(this.noZoom || this.state != this.STATE.TOUCH_DOLLY) return;
            let dx = event.touches[0].clientX - event.touches[1].clientX;
            let dy = event.touches[0].clientY - event.touches[1].clientY;
            let distance = Math.sqrt(dx * dx + dy * dy);
            this.dollyEnd.set(0, distance);
            this.dollyDelta.subVectors(this.dollyEnd, this.dollyStart);
            if(this.dollyDelta.y > 0) this.dollyOut();
            else this.dollyIn();
            this.dollyStart.copy(this.dollyEnd);
            break;
        case 3:
            // three-fingered touch: pan
            if(this.noPan || this.state != this.STATE.TOUCH_PAN) return;
            this.panEnd.set(event.touches[0].clientX, event.touches[0].clientY);
            this.panDelta.subVectors(this.panEnd, this.panStart);
            this.pan(this.panDelta.x, this.panDelta.y);
            this.panStart.copy(this.panEnd);
            break;
        default:
            this.touchEnd();
            return;
        }

        this.update();
    }

    touchEnd() {
        if(!this.enabled) return;
        let de = document.documentElement;
        de.removeEventListener('touchend',    this.touchEnd,  false);
        de.removeEventListener('touchmove',   this.touchMove, false);
        de.removeEventListener('touchcancel', this.touchEnd,  false);
        this.dispatchEvent(this.endEvent);
        this.state = this.STATE.NONE;
    }
}



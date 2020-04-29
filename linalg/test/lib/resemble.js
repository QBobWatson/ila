
import should from 'should';
import Matrix from '../../src/matrix.js';


should.use(function(should, Assertion) {
    // Custom assertion for comparing arrays of numbers using 'approximately'
    function approximate(obj, other, ε) {
        if(typeof obj === "number")
            obj.should.be.approximately(other, ε);
        else if(obj instanceof Matrix) {
            other.should.be.an.instanceOf(Matrix);
            obj.m.should.equal(other.m);
            obj.n.should.equal(other.n);
            approximate([...obj.rows()], [...other.rows()], ε);
        }
        else {
            obj.should.be.an.Array();
            other.should.be.an.Array();
            other.should.have.length(obj.length);
            for(let i = 0; i < obj.length; ++i)
                approximate(obj[i], other[i], ε);
        }
    }

    Assertion.add('resemble', function(other, ε=1e-10) {
        this.params = {
            operator: 'to approximately equal',
            expected: other
        };
        approximate(this.obj, other, ε);
    });
});

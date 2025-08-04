package abyss.inference;

import beast.base.evolution.likelihood.TreeLikelihood;

public class CachedTreeLikelihood extends TreeLikelihood {

    protected double[] storedPatternLogLikelihoods;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        patternLogLikelihoods = new double[alignment.getPatternCount()];
        storedPatternLogLikelihoods = new double[alignment.getPatternCount()];
    }

    
    @Override
    public double[] getPatternLogLikelihoods() {
        if (beagle != null && isDirtyCalculation()) {
            System.arraycopy(beagle.getPatternLogLikelihoods(), 0, patternLogLikelihoods, 0, patternLogLikelihoods.length);
        }
		return patternLogLikelihoods.clone();
    }
    
    @Override
    public void store() {
        super.store();
        System.arraycopy(patternLogLikelihoods, 0, storedPatternLogLikelihoods, 0, patternLogLikelihoods.length);
    }

    @Override
    public void restore() {
        super.restore();
        double[] tmp = patternLogLikelihoods;
        patternLogLikelihoods = storedPatternLogLikelihoods;
        storedPatternLogLikelihoods = tmp;
    }

}

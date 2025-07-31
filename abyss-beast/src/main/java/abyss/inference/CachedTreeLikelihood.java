package abyss.inference;

import beast.base.evolution.likelihood.TreeLikelihood;

public class CachedTreeLikelihood extends TreeLikelihood {

    protected double[] storedPatternLogLikelihoods;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        storedPatternLogLikelihoods = new double[alignment.getPatternCount()];
    }

    @Override
    public void store() {
        super.store();
        if (beagle != null) return;
        System.arraycopy(patternLogLikelihoods, 0, storedPatternLogLikelihoods, 0, patternLogLikelihoods.length);
    }

    @Override
    public void restore() {
        super.restore();
        if (beagle != null) return;
        double[] tmp = patternLogLikelihoods;
        patternLogLikelihoods = storedPatternLogLikelihoods;
        storedPatternLogLikelihoods = tmp;
    }

}

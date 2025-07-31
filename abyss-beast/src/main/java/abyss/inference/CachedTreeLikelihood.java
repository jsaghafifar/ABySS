package abyss.inference;

import beast.base.evolution.likelihood.TreeLikelihood;
import beast.pkgmgmt.BEASTClassLoader;

import java.lang.reflect.InvocationTargetException;

public class CachedTreeLikelihood extends TreeLikelihood {

    protected double[] storedPatternLogLikelihoods;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        try {
            Object o = newTreeLikelihood();
            if (o instanceof CachedBeagleTreeLikelihood) {
                beagle = (CachedBeagleTreeLikelihood) o;
            }
            beagle.initByName(
                    "data", dataInput.get(), "tree", treeInput.get(), "siteModel", siteModelInput.get(),
                    "branchRateModel", branchRateModelInput.get(), "useAmbiguities", m_useAmbiguities.get(),
                    "useTipLikelihoods", m_useTipLikelihoods.get(),"scaling", scaling.get().toString(),
                    "rootFrequencies", rootFrequenciesInput.get());
            if (beagle.getBeagle() != null) {
                //a Beagle instance was found, so we use it
                return;
            }
        } catch (Exception e) {
            // ignore
        }
        // No Beagle instance was found, so we use java
        beagle = null;

        storedPatternLogLikelihoods = new double[alignment.getPatternCount()];
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

    private TreeLikelihood newTreeLikelihood() throws InstantiationException, IllegalAccessException, IllegalArgumentException, InvocationTargetException, NoSuchMethodException, SecurityException, ClassNotFoundException {
        // must set implementationInput to CachedBeagleTreeLikelihood to use beagle
        String class_ = implementationInput.get();
        return (TreeLikelihood) BEASTClassLoader.forName(class_).getConstructor().newInstance();
    }

}

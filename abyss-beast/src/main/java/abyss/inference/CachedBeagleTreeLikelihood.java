package abyss.inference;

import beagle.Beagle;
import beast.base.evolution.likelihood.BeagleTreeLikelihood;

/**
 * @author Remco Bouckaert
 */
public class CachedBeagleTreeLikelihood extends BeagleTreeLikelihood {

    @Override
    public void restore() {
        super.restore();
        
        double[] sumLogLikelihoods = new double[1];
        int rootIndex = partialBufferHelper.getOffsetIndex(treeInput.get().getRoot().getNr());

        int cumulateScaleBufferIndex = Beagle.NONE;
        if (useScaleFactors) {
           cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
        }
        
        getBeagle().calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0},
                new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);
    }

}

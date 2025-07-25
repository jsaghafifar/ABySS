package abyss.substitutionmodel;

import abyss.distributions.MixedTreeLikelihood;
import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.LikelihoodCore;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import org.apache.commons.math3.util.FastMath;

import java.io.PrintStream;
import java.util.Arrays;

@Description("Site logger for models with a number of mixed likelihoods.")
public class SiteModelLogger extends BEASTObject implements Loggable {

    final public Input<MixedTreeLikelihood> mixedLikelihoodsInput = new Input<>("mixedLikelihoods",
            "Mixed tree likelihoods", Input.Validate.REQUIRED);

    int siteCount;
    int nrOfPatterns;
    int modelCount;

    @Override
    public void initAndValidate() {
        this.modelCount = mixedLikelihoodsInput.get().pLikelihoods.get().size();
        Alignment data = (Alignment) mixedLikelihoodsInput.get().pLikelihoods.get().get(0).getInput("data").get();
        this.siteCount = data.getSiteCount();
        this.nrOfPatterns = data.getPatternCount();
    }


    @Override
    public void init(PrintStream out) {
        for (int i = 0; i < siteCount; i++) {
            out.print((getID() != null ? getID() : "site") + "." + (i+1) + "\t");
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        // sample from models posteriors
        double[] partialLogLikelihoods = mixedLikelihoodsInput.get().getPartialLogLikelihoods();
        double[] posteriorOfEachModel = new double[partialLogLikelihoods.length];
        double max = Arrays.stream(partialLogLikelihoods).max().getAsDouble();
        double pSum = 0;
        for (int i = 0; i < partialLogLikelihoods.length; i++) {
            posteriorOfEachModel[i] = FastMath.exp(partialLogLikelihoods[i]-max);
            pSum += posteriorOfEachModel[i];
        }
        for (int i = 0; i < posteriorOfEachModel.length; i++) {
            posteriorOfEachModel[i] /= pSum;
        }
        int modelIndex = Randomizer.randomChoicePDF(posteriorOfEachModel);

        String mode = mixedLikelihoodsInput.get().modeInput.get();
        if (!mode.equalsIgnoreCase("mix") && modelIndex < modelCount) {
            for (int i = 0; i < siteCount; i++) {
                out.print(modelIndex + "\t");
            }
        } else {
            // for site mixture, calculate likelihood for every site, sampling model index from each
            double[] siteModelWeights = mixedLikelihoodsInput.get().siteModelWeightsInput.get().getDoubleValues();
            double[][] rootPartials = new double[this.nrOfPatterns][this.modelCount];

            // get root partials for site likelihoods
            for (int i = 0; i < modelCount; i++) {
                TreeLikelihood likelihood = (TreeLikelihood) mixedLikelihoodsInput.get().pLikelihoods.get().get(i);
                LikelihoodCore core = likelihood.getLikelihoodCore();
                Tree tree = (Tree) likelihood.treeInput.get();
                core.getNodePartials(tree.getRoot().getNr(), rootPartials[i]);
            }

            double[][] posteriorOfEachModelPerPattern = new double[this.nrOfPatterns][this.modelCount];
            double[] probsPattern = new double[this.modelCount];

            for (int i = 0; i < this.nrOfPatterns; i++) {
                pSum = 0;
                for (int j = 0; j < this.modelCount; j++) {
                    probsPattern[j] = rootPartials[i][j] * siteModelWeights[j];
                    pSum += probsPattern[j];
                }

                // normalise posterior
                for (int j = 0; j < this.modelCount; j++) {
                    posteriorOfEachModelPerPattern[i][j] = probsPattern[j] / pSum;
                }
            }

            Alignment data = ((TreeLikelihood) mixedLikelihoodsInput.get().pLikelihoods.get().get(0)).dataInput.get();

            for (int i = 0; i < siteCount; i++) {
                int patternNum = data.getPatternIndex(i);
                double[] p = posteriorOfEachModelPerPattern[patternNum];
                double sum = 0;
                for (int j = 0; j < p.length; j++) {
                    sum += p[j];
                }
                if (Math.abs(sum-1) > 1e-6) out.print(-1 + "\t");
                else out.print(Randomizer.randomChoicePDF(p) + "\t");
            }
        }
    }

    @Override
    public void close(PrintStream out) {

    }
}

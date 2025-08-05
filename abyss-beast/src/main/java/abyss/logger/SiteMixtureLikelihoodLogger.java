package abyss.logger;

import abyss.distributions.MixedTreeLikelihood;
import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;

import java.io.PrintStream;

@Description("Model logger for models with a number of mixed likelihoods.")
public class SiteMixtureLikelihoodLogger extends BEASTObject implements Loggable {

    final public Input<MixedTreeLikelihood> mixedLikelihoodsInput = new Input<>("mixedLikelihoods",
            "Mixed tree likelihoods", Input.Validate.REQUIRED);

    int siteMixtureIndex;

    @Override
    public void initAndValidate() {
        if (mixedLikelihoodsInput.get().modeInput.equals("avg"))
            throw new IllegalArgumentException("No site mixture models involved in 'avg' mode.");
        this.siteMixtureIndex = mixedLikelihoodsInput.get().modeInput.equals("both") ?
                mixedLikelihoodsInput.get().pLikelihoods.get().size() : 0;
    }


    @Override
    public void init(PrintStream out) {
        out.print((getID() != null ? getID() : "siteMixtureLikelihood") + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        // get site mixture log likelihood
        double likelihood = mixedLikelihoodsInput.get().getPartialLogLikelihoods()[this.siteMixtureIndex];
        // if model averaging involved, need to remove meta weight
        if (mixedLikelihoodsInput.get().modeInput.equals("both"))
            likelihood -= Math.log(mixedLikelihoodsInput.get().getMetaWeights()[this.siteMixtureIndex]);
        out.print(likelihood + "\t");

    }

    @Override
    public void close(PrintStream out) {

    }
}

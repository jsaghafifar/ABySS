package abyss.lphybeast.tobeast.generator;

import abyss.MixedAlignment;
import abyss.distributions.MixedTreeLikelihood;
import abyss.inference.CachedTreeLikelihood;
import abyss.logger.AlignmentModelLogger;
import abyss.logger.SiteMixtureLikelihoodLogger;
import abyss.logger.SiteModelLogger;
import beast.base.core.BEASTInterface;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.inference.distribution.Dirichlet;
import beast.base.inference.distribution.Prior;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.likelihood.PhyloCTMC;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import java.util.ArrayList;
import java.util.List;

import static lphybeast.BEASTContext.createPrior;
import static lphybeast.tobeast.generators.PhyloCTMCToBEAST.constructSiteModel;
import static lphybeast.tobeast.generators.PhyloCTMCToBEAST.constructTreeAndBranchRate;

/**
 * @author Jasmine Saghafifar
 */
public class MixedAlignmentToBEAST implements GeneratorToBEAST<MixedAlignment, MixedTreeLikelihood> {

    public MixedTreeLikelihood generatorToBEAST(MixedAlignment mixedAlignment, BEASTInterface value, BEASTContext context) {
        return createMixedTreeLikelihood(mixedAlignment, value, context);
    }

    private MixedTreeLikelihood createMixedTreeLikelihood(MixedAlignment mixedAlignment, BEASTInterface value, BEASTContext context) {
        assert value instanceof beast.base.evolution.alignment.Alignment;
        beast.base.evolution.alignment.Alignment alignment = (beast.base.evolution.alignment.Alignment)value;

        MixedTreeLikelihood treeLikelihood = new MixedTreeLikelihood();

        Value<Alignment> aln1 = mixedAlignment.getAlignment1();
        Value<Alignment> aln2 = mixedAlignment.getAlignment2();

        // TODO make sure aln1 and aln2 and any other alignments are not included in XML

        List<PhyloCTMC> phyloCTMCList = new ArrayList<>();
        phyloCTMCList.add((PhyloCTMC) aln1.getGenerator());
        phyloCTMCList.add((PhyloCTMC) aln2.getGenerator());

        List<CachedTreeLikelihood> likelihoods = new ArrayList<>();
        int i=0;
        // initialise likelihoods within
        for (PhyloCTMC phyloCTMC : phyloCTMCList) {
            CachedTreeLikelihood likelihood = new CachedTreeLikelihood();

            constructTreeAndBranchRate(phyloCTMC, likelihood, context);
            SiteModel siteModel = constructSiteModel(phyloCTMC, context);

            likelihood.setInputValue("siteModel", siteModel);
            likelihood.setInputValue("data", alignment);
            likelihood.initAndValidate();

            likelihood.setID(alignment.getID() + ".treeLikelihood"+i);
            context.addExtraLoggable(likelihood);
            likelihoods.add(likelihood);
            i++;
        }

        treeLikelihood.setInputValue("likelihood", likelihoods);
        treeLikelihood.setInputValue("mode", "both");

        addSiteMixtureWeightsPrior(context,mixedAlignment,treeLikelihood);

        treeLikelihood.initAndValidate();
        treeLikelihood.setID(alignment.getID() + ".treeLikelihood");

//        not needed?
        context.addExtraLoggable(treeLikelihood);

        // loggers
        addSiteMixtureLikelihoodLogger(context,treeLikelihood);
        addModelIndicatorLogger(context,treeLikelihood);
        addSiteLikelihoodsLogger(context,treeLikelihood);

        return treeLikelihood;
    }

    private void addSiteMixtureWeightsPrior(BEASTContext context, MixedAlignment alignment, MixedTreeLikelihood treeLikelihood) {
        Value<Double[]> weights;

        if (alignment.getSiteMixtureWeights() != null)
            weights = alignment.getSiteMixtureWeights();
        else throw new IllegalArgumentException("Site mixture weights must be specified for lphybeast.");

        beast.base.inference.distribution.Dirichlet dirichlet = new Dirichlet();
        dirichlet.setInputValue("alpha", context.getAsRealParameter(weights));

        Prior prior = createPrior(dirichlet, context.getAsRealParameter(weights));
        context.addBEASTObject(prior, weights);
        context.addExtraLoggable(prior);
        treeLikelihood.setInputValue("siteModelWeights", context.getAsRealParameter(weights));
    }

    private void addSiteMixtureLikelihoodLogger(BEASTContext context, MixedTreeLikelihood treeLikelihood) {
        SiteMixtureLikelihoodLogger mixtureLikelihoodLogger = new SiteMixtureLikelihoodLogger();
        mixtureLikelihoodLogger.setInputValue("mixedLikelihoods", treeLikelihood);
        mixtureLikelihoodLogger.initAndValidate();
        context.addExtraLoggable(mixtureLikelihoodLogger);
    }

    private void addModelIndicatorLogger(BEASTContext context, MixedTreeLikelihood treeLikelihood) {
        AlignmentModelLogger indicatorLogger = new AlignmentModelLogger();
        indicatorLogger.setInputValue("mixedLikelihoods", treeLikelihood);
        indicatorLogger.initAndValidate();
        context.addExtraLoggable(indicatorLogger);
    }

    private void addSiteLikelihoodsLogger(BEASTContext context, MixedTreeLikelihood treeLikelihood) {
        SiteModelLogger siteLogger = new SiteModelLogger();
        siteLogger.setInputValue("mixedLikelihoods", treeLikelihood);
        siteLogger.initAndValidate();
        context.addExtraLoggable(siteLogger);
    }

    @Override
    public Class<MixedAlignment> getGeneratorClass() {
        return MixedAlignment.class;
    }

    @Override
    public Class<MixedTreeLikelihood> getBEASTClass() {
        return MixedTreeLikelihood.class;
    }
}
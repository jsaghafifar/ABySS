open module abyss.beast {
    requires beast.pkgmgmt;
    requires beast.base;

    requires beagle;
//    requires beast.fx;
    requires colt;
//    requires commons.math3;
    requires beast.classic;
    requires org.apache.commons.math4.core;
    requires org.apache.commons.math4.legacy;
    requires org.apache.commons.numbers.gamma;
    requires org.apache.commons.statistics.distribution;
    requires jdk.jfr;

    exports abyss.distributions;
    exports abyss.logger;
    exports abyss.inference;
    exports abyss.substitutionmodel;

    provides beast.base.core.BEASTInterface with
            abyss.distributions.InformedDirichletPrior,
            abyss.distributions.MixedTreeLikelihood,
            abyss.distributions.SVSPrior,
            abyss.distributions.PseudoPrior,
            abyss.distributions.EigenFriendlyQPrior,
            abyss.logger.ABySSFrequencyLogger,
            abyss.logger.RootMeanSquareLogger,
            abyss.logger.AlignmentModelLogger,
            abyss.logger.NetFluxLogger,
            abyss.logger.DetailedBalanceLogger,
            abyss.logger.SiteMixtureLikelihoodLogger,
            abyss.logger.SiteModelLogger,
            abyss.logger.ClassicDetailedBalanceLogger,
            abyss.logger.ClassicFrequencyLogger,
            abyss.logger.ClassicNetFluxLogger,
            abyss.inference.CachedTreeLikelihood,
            abyss.inference.CachedBeagleTreeLikelihood,
            abyss.substitutionmodel.ABySSubstitutionModel,
            abyss.substitutionmodel.ABySSModelAveraging;
}
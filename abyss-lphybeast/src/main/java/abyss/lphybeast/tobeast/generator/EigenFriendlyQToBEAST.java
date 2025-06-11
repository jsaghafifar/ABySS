package abyss.lphybeast.tobeast.generator;

import abyss.NonReversible;
import abyss.distributions.EigenFriendlyQPrior;
import beast.base.core.BEASTInterface;
import beast.base.inference.parameter.RealParameter;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

public class EigenFriendlyQToBEAST implements GeneratorToBEAST<NonReversible, EigenFriendlyQPrior> {
    @Override
    public EigenFriendlyQPrior generatorToBEAST(NonReversible generator, BEASTInterface value, BEASTContext context) {
        if (generator.getSymmetric().value()) return null;

        EigenFriendlyQPrior qPrior = new EigenFriendlyQPrior();

        Value<Double[]> rates = generator.getRates();
        RealParameter ratesParameter = (RealParameter)context.getBEASTObject(rates);
        qPrior.setInputValue("rates", ratesParameter);

        if (generator.getIndicators() != null) {
            RealParameter rateIndicatorParameter = (RealParameter)context.getBEASTObject(generator.getIndicators());
            qPrior.setInputValue("rateIndicator", rateIndicatorParameter);
        }

        double root = (-1 - Math.sqrt(1+4* rates.value().length))/2;
        int numStates = (int) Math.abs(root);
        qPrior.setInputValue("nrOfStates", numStates);

        qPrior.initAndValidate();
        return qPrior; // TODO fix so lphybeast doesn't treat as substModel to init
    }

    @Override
    public Class<NonReversible> getGeneratorClass() {
        return NonReversible.class;
    }

    @Override
    public Class<EigenFriendlyQPrior> getBEASTClass() {
        return EigenFriendlyQPrior.class;
    }
}

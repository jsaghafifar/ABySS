package abyss.distributions;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.distribution.Prior;
import beast.base.inference.parameter.IntegerParameter;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.Random;

@Description("Pseudo prior to be used when a model is not active in model averaging.")
public class PseudoPrior extends Distribution {

    final public Input<Prior> priorInput = new Input<>("prior", "distribution to be used when model is on.", Input.Validate.REQUIRED);
    final public Input<Prior> pseudoPriorInput = new Input<>("pseudo", "informed distribution to be used when model is off.", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> modelIndicatorInput = new Input<>("modelIndicator", "the active model.", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> modelIndexInput = new Input<>("modelIndex", "the prior's model(s), as its index in the model list.", Input.Validate.REQUIRED);

    public double calculateLogP() {
        this.logP = 0.0;
        boolean priorActive = false;
        if (modelIndexInput.get().getDimension() == 1) {
            priorActive = Objects.equals(modelIndexInput.get().getValue(), modelIndicatorInput.get().getValue());
        } else {
            for (int i = 0; i < modelIndexInput.get().getDimension(); i++) {
                if (modelIndicatorInput.get().getValue() == modelIndexInput.get().getValue(i)) {
                    priorActive = true;
                    break;
                }
            }
        }

        Prior prior;
        if (priorActive) prior = priorInput.get();
        else prior = pseudoPriorInput.get();

        logP = prior.calculateLogP();

        return logP;
    }

    @Override
    public List<String> getArguments() {
        List<String> args = new ArrayList<>();
        args.add(modelIndexInput.get().getID());
        return args;
    }

    @Override
    public List<String> getConditions() {
        List<String> conds = new ArrayList<>();
        conds.add(priorInput.get().getID());
        conds.add(pseudoPriorInput.get().getID());
        conds.add(modelIndicatorInput.get().getID());
        return conds;
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException();
    }

}

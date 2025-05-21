package abyss;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;

//TODO flu model citation
public class FLURates<U> extends DeterministicFunction<U[]> {

    public FLURates() {}

    @GeneratorInfo(name = "fluRates", description = "")
    public Value<U[]> apply() {
        FLU flu = new FLU(null);
        U[] array = (U[]) flu.getRates();
        return new Value<>( array, this);
    }
}
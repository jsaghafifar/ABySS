package abyss;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;

//TODO flu model citation
public class FLUFreqs<U> extends DeterministicFunction<U[]> {

    public FLUFreqs() {}

    @GeneratorInfo(name = "fluFreqs", description = "")
    public Value<U[]> apply() {
        FLU flu = new FLU(null);
        U[] array = (U[]) flu.getFreq();
        return new Value<>( array, this);
    }
}
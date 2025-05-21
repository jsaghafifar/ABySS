package abyss;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;

//TODO model citation
public class HIVBRates<U> extends DeterministicFunction<U[]> {

    public HIVBRates() {}

    @GeneratorInfo(name = "hivbRates", description = "")
    public Value<U[]> apply() {
        HIVB hivb = new HIVB(null);
        U[] array = (U[]) hivb.getRates();
        return new Value<>( array, this);
    }
}
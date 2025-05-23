package abyss;

import lphy.base.evolution.substitutionmodel.RateMatrix;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.model.datatype.DoubleArray2DValue;

//TODO citation here
public class HIVB extends RateMatrix {

    public HIVB(@ParameterInfo(name = RateMatrix.meanRateParamName, description = "the mean rate of the process. default 1.0", optional=true) Value<Number> meanRate) {
        super(meanRate);
    }

    @GeneratorInfo(name = "hivb", verbClause = "is", narrativeName = "Empirical HIVb amino acid substitution model",
            category = GeneratorCategory.RATE_MATRIX,
            description = "A symmetric instantaneous rate matrix for HIVb amino acid substitution.")
    public Value<Double[][]> apply() {
        Double[][] Q = getQ();
        return new DoubleArray2DValue(Q, this);
    }

//    @Override
//    public boolean canReturnComplexDiagonalization() { return true; }

    protected Double[][] getQ() {
        Double[] r = getRates();
        Double[] f = getFreq();
        int numStates = f.length;
        Double[][] Q = new Double[numStates][numStates];
        int x = 0;
        for (int i = 0; i < f.length; i++) {
            for (int j = i + 1; j < f.length; j++) {
                Q[i][j] = r[x] * f[j];
                Q[j][i] = r[x] * f[i];
                x++;
            }
            Q[i][i] = 0.0;
            for (int j = 0; j < f.length; j++) {
                if (i==j) continue;
                Q[i][i] -= Q[i][j];
            }
        }
        normalize(f, Q, totalRateDefault1());
        return Q;
    }

    public Double[] getRates() {
        Double[] rates;
        rates = new Double[]{0.123758000,
                             1.455040000,0.005000000,
                             1.481350000,0.005000000,10.58720000,
                             0.014126900,9.298150000,0.005000000,0.005000000,
                             2.135360000,0.897871000,2.838060000,3.927750000,0.291561000,
                             0.084761300,0.240073000,1.916900000,0.119740000,0.145558000,0.005000000,
                             0.005000000,0.005000000,0.017679200,0.006090790,3.398360000,0.005000000,0.103111000,
                             0.005000000,0.005000000,0.005000000,4.614820000,0.034265800,0.521705000,0.005000000,0.322319000,
                             0.215256000,0.129777000,0.008760480,0.005000000,8.524840000,0.005000000,1.741710000,5.958790000,0.081499500,
                             0.018664300,0.005000000,0.005000000,0.175789000,0.188025000,0.005000000,0.005000000,11.20650000,1.282460000,5.319610000,
                             0.005000000,0.086064200,17.66120000,0.079263300,0.005000000,0.323401000,7.645850000,0.680565000,7.904430000,0.005000000,0.005000000,
                             2.122170000,0.005000000,0.034265800,0.012022600,0.005000000,0.005000000,2.453180000,0.041059300,0.031386200,2.077570000,0.005000000,0.007395780,
                             0.055112800,0.005000000,0.005000000,2.560200000,0.005000000,0.061913700,7.055450000,0.005000000,6.547370000,1.494560000,0.303676000,0.672052000,4.472110000,
                             0.307507000,0.351721000,0.005000000,0.074921800,0.005000000,3.653450000,9.040440000,0.677289000,20.45000000,0.701427000,2.513940000,0.295543000,1.283550000,3.421500000,
                             2.466330000,4.693140000,0.528230000,0.005000000,0.956472000,4.380410000,0.382747000,1.218030000,0.504111000,0.927656000,0.005000000,13.14470000,5.377620000,0.116311000,3.479100000,
                             15.91830000,0.739969000,0.274724000,0.289774000,0.014126900,0.369615000,0.711594000,8.612170000,4.671420000,0.043767300,4.940260000,6.886670000,2.014170000,0.243589000,2.868680000,8.931070000,
                             7.614280000,0.420027000,1.047930000,1.028470000,0.723274000,0.953155000,0.005000000,17.73890000,0.265829000,1.410360000,6.853200000,0.026656000,0.005000000,0.020915300,0.081245400,0.074921800,0.709226000,
                             0.005000000,2.632770000,0.005000000,0.005000000,0.829343000,1.216740000,0.069517900,0.005000000,0.005000000,0.748843000,0.089078000,0.005000000,0.044450600,0.026656000,0.991338000,0.024872800,0.005000000,0.005000000,
                             0.005000000,7.579320000,0.674653000,0.079263300,15.34000000,0.005000000,18.69430000,0.148168000,0.005000000,0.111986000,0.005000000,1.764170000,0.030438100,0.113033000,0.009918260,0.648024000,0.105652000,0.041059300,1.280220000};
        double sum = 0;
        for (int i = 0; i < rates.length; i++) {
            sum += rates[i];
        }
        for (int i = 0; i < rates.length; i++) {
            rates[i] /= sum;
        }
        return rates;
    }

    public Double[] getFreq() {
        Double[] freq = new Double[20]; // ORDER: ACDEFGHIKLMNPQRSTVWY
        freq[0] = 0.060490222;
        freq[1] = 0.020075899;
        freq[2] = 0.042109048;
        freq[3] = 0.071567447;
        freq[4] = 0.028809447;
        freq[5] = 0.072308239;
        freq[6] = 0.022293943;
        freq[7] = 0.069730629;
        freq[8] = 0.056968211;
        freq[9] = 0.098851122;
        freq[10]= 0.019768318;
        freq[11]= 0.044127815;
        freq[12]= 0.046025282;
        freq[13]= 0.053606488;
        freq[14]= 0.066039665;
        freq[15]= 0.050604330;
        freq[16]= 0.053636813;
        freq[17]= 0.061625237;
        freq[18]= 0.033011601;
        freq[19]= 0.028350243;
        double sum = 0;
        for (int i = 0; i < freq.length; i++) {
            sum += freq[i];
        }
        for (int i = 0; i < freq.length; i++) {
            freq[i] /= sum;
        }
        return freq;
    }

}
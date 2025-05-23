package abyss;

import lphy.base.evolution.substitutionmodel.RateMatrix;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.model.datatype.DoubleArray2DValue;

//TODO citation here
public class HIVW extends RateMatrix {

    public HIVW(@ParameterInfo(name = RateMatrix.meanRateParamName, description = "the mean rate of the process. default 1.0", optional=true) Value<Number> meanRate) {
        super(meanRate);
    }

    @GeneratorInfo(name = "hivw", verbClause = "is", narrativeName = "Empirical HIVw amino acid substitution model",
            category = GeneratorCategory.RATE_MATRIX,
            description = "A symmetric instantaneous rate matrix for HIVw amino acid substitution.")
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
        rates = new Double[]{0.167653000,
                             4.435210000,0.005000000,
                             5.563250000,0.005000000,12.12330000,
                             0.597923000,0.362959000,0.005000000,0.005000000,
                             1.868500000,0.048979800,10.39690000,14.78010000,0.005000000,
                             0.005000000,0.005000000,2.317790000,0.005000000,0.005000000,0.005000000,
                             0.005000000,0.005000000,0.145124000,0.039051200,1.482880000,0.005000000,0.005000000,
                             0.592784000,0.005000000,0.894313000,23.96260000,0.005000000,0.279425000,0.224060000,0.817481000,
                             0.160240000,0.005000000,0.005000000,0.129839000,7.487810000,0.048979800,1.763820000,9.102460000,0.005000000,
                             0.005000000,0.005000000,0.005000000,0.005000000,0.005000000,0.048979800,0.005000000,17.30640000,4.095640000,11.38390000,
                             0.617509000,0.060493200,29.40870000,0.201526000,0.005000000,0.060493200,8.598760000,0.987028000,10.66550000,0.005000000,0.201526000,
                             1.009810000,0.005000000,0.005000000,0.005000000,0.034225200,0.005000000,13.94440000,0.005000000,0.111928000,9.830950000,0.005000000,0.344848000,
                             0.005000000,0.005000000,0.005000000,3.206560000,0.005000000,0.060493200,18.54650000,0.034225200,13.07050000,2.890480000,0.005000000,0.342068000,3.045020000,
                             0.074480800,2.863640000,0.067453900,0.025163200,0.005000000,13.43790000,6.844050000,1.340690000,39.88970000,0.586757000,3.286520000,0.160240000,0.404723000,10.67460000,
                             8.594200000,1.121950000,0.427881000,0.005000000,4.279390000,6.279660000,0.725157000,0.740091000,0.005000000,6.143960000,0.392575000,14.56990000,14.24900000,0.160240000,8.350240000,
                             24.14220000,0.005000000,0.630395000,0.458743000,0.114512000,0.048979800,0.959560000,9.363450000,4.048020000,0.005000000,7.413130000,4.542060000,4.337010000,0.203091000,0.928203000,6.340790000,
                             24.80940000,0.005000000,2.917860000,2.199520000,2.280000000,2.796220000,0.827479000,24.82310000,0.128065000,2.953440000,14.76830000,0.074480800,0.005000000,0.005000000,0.279425000,0.862637000,0.005000000,
                             0.005000000,5.498940000,0.005000000,0.005000000,0.005000000,2.825800000,0.005000000,0.005000000,0.005000000,1.370310000,0.005000000,0.005000000,0.005000000,0.044329800,5.965640000,1.101560000,0.005000000,0.005000000,
                             0.005000000,8.348350000,2.281540000,0.005000000,4.127280000,0.005000000,47.48890000,0.114512000,0.005000000,0.005000000,0.579198000,5.064750000,0.005000000,0.005000000,0.005000000,0.933142000,0.490608000,1.354820000,0.005000000};
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
        freq[0] = 0.0377494;
        freq[1] = 0.0240105;
        freq[2] = 0.0342034;
        freq[3] = 0.0618606;
        freq[4] = 0.0422741;
        freq[5] = 0.0838496;
        freq[6] = 0.0156076;
        freq[7] = 0.0983641;
        freq[8] = 0.0641682;
        freq[9] = 0.0577867;
        freq[10]= 0.0158419;
        freq[11]= 0.0891129;
        freq[12]= 0.0458601;
        freq[13]= 0.0437824;
        freq[14]= 0.0573210;
        freq[15]= 0.0550846;
        freq[16]= 0.0813774;
        freq[17]= 0.0515638;
        freq[18]= 0.0195970;
        freq[19]= 0.0205847;
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
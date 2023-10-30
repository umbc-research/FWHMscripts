from scipy.optimize import curve_fit
from numpy import exp, linspace, random, indices,bincount, sqrt

def gaussian_1d(x, mu, sigma, amplitude, offset):
    return amplitude * exp( -((x-mu)/sigma)**2/2 ) + offset


def fit_gaussian_1d(x_data, y_data, bounds=([40, 1, 500, 1],[60,100,75000, 25000])):
    #bounds go like this: ([minMu, minSig, minAmp, minOff],[maxMu, maxSig, maxAmp, maxOff])
    params, cov = curve_fit(gaussian_1d, x_data, y_data, bounds=([40, 1, 500, 1],[60,100,75000, 25000]))
    return params

def extract_radial_data(data, xC, yC):
    #Get matrix of integer indices associated with subFrame
    y, x = indices((data.shape))

    #Generate matrix of radius values
    r = sqrt((x - xC)**2 + (y - yC)**2)

    #Force integer values (np.sqrt gives floats)
    r = r.astype(int)

    #Generate a histogram of radius bin, each weighed by corresponding counts
    weightedRadiusHistogram = bincount(r.ravel(), weights=data.ravel())
    unweightedRadiusHistogram = bincount(r.ravel())

    #Get average for each radius bin
    averageCountsPerRadiusBin = weightedRadiusHistogram / unweightedRadiusHistogram
    return averageCountsPerRadiusBin




#code for testing functions in profileFitting

if __name__ == "__main__":
    from sys import argv
    from matplotlib import pyplot as plt
    
    mu, sigma, amplitude, offset = float(argv[1]), float(argv[2]), float(argv[3]), float(argv[4]) 
    x = linspace(0,100,1000)
    y = gaussian_1d(x, mu, sigma, amplitude, offset)

    y_rand = y + random.rand(len(y))*2.5 - 1.25
    m, s, a, o = fit_gaussian_1d(x[::25], y_rand[::25])
    plt.plot(x[::25], y_rand[::25], 'r.')
    plt.plot(x[::25], gaussian_1d(x[::25], m, s, a, o), 'b')
    print(m, s, a, o)
    plt.grid(1)
    plt.show()
            
    
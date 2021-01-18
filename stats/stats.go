package stats

import (
    "fmt"
    "math"
    . "github.com/jhpeng/ctQMC/dtype"
)

func New() Estimator {
    var est Estimator

    est.Sample = make([]float64,100000)
    est.N = 0

    return est
}

func Append(est Estimator, sample float64) Estimator {
    if (est.N+1)>len(est.Sample) {
        temp := make([]float64,len(est.Sample)*2)
        copy(temp,est.Sample)
        est.Sample = temp
    }

    est.Sample[est.N] = sample
    est.N++

    return est
}

func Cut(est Estimator, length int) Estimator {
    temp := make([]float64,len(est.Sample))

    if length<est.N {
        n:=0
        for i:=length;i<est.N;i++ {
            temp[n] = est.Sample[i]
            n++
        }
        est.N = n
    } else {
        est.N = 0
    }
    est.Sample = temp

    return est
}

func Mean(est Estimator) float64 {
    var mean float64 = 0

    for i:=0;i<est.N;i++ {
        mean += est.Sample[i]
    }
    mean = mean/float64(est.N)

    return mean
}

func Std(est Estimator) float64 {
    var std float64 = 0
    var div float64 = 0
    mean := Mean(est)

    for i:=0;i<est.N;i++ {
        div = est.Sample[i]-mean
        std += div*div
    }
    std = math.Sqrt(std/float64(est.N))

    return std
}

// <x_i * x_{i+dist}>
// n = N-dist
// 1/n sum_{i=0}^n x_i * x_{i+dist}
func Correlation(est Estimator, dist int) float64 {
    var cross float64 = 0
    n := est.N-dist

    if n<=0 {
        fmt.Printf("Estimator.N must be larger than dist!\n")
    }

    for i:=0;i<n;i++ {
        cross += est.Sample[i]*est.Sample[i+dist]
    }
    cross = cross/float64(n)

    return cross
}

func AutocorrelationTime(est Estimator) float64 {
    var tau float64 = 0

    if est.N<1000 {
        fmt.Printf("Estimator.N must be larger than 1000!\n")
    }

    mean := Mean(est)

    for i:=1;i<10;i++ {
        chi_k := Correlation(est,i) - mean*mean
        tau += (1.0-float64(i)/float64(est.N))*chi_k
    }
    tau = tau/(Correlation(est,0)-mean*mean)

    return tau
}

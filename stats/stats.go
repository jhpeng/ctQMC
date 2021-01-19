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
    est.Block  = make([]float64,1024)
    est.Bsize  = 2
    est.Nblock = 0

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

    if est.N==est.Bsize*(est.Nblock+1) {
        var b float64 = 0
        for i:=0;i<est.Bsize;i++ {
            b += est.Sample[est.Nblock*est.Bsize+i]
        }
        b = b/float64(est.Bsize)
        est.Block[est.Nblock] = b
        est.Nblock++

        if est.Nblock==1024 {
            for i:=0;i<512;i++ {
                est.Block[i] = est.Block[2*i]+est.Block[2*i+1]
                est.Block[i] = est.Block[i]*0.5
            }
            est.Nblock = 512
            est.Bsize *= 2
        }
    }

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

func MeanBlock(est Estimator) float64 {
    var mean float64 = 0

    for i:=0;i<est.Nblock;i++ {
        mean += est.Block[i]
    }
    mean = mean/float64(est.Nblock)

    return mean
}

func StdBlock(est Estimator) float64 {
    var std float64 = 0
    var div float64 = 0
    mean := MeanBlock(est)

    for i:=0;i<est.Nblock;i++ {
        div = est.Block[i]-mean
        std += div*div
    }
    std = math.Sqrt(std/float64(est.Nblock))

    return std
}

func StdError(est Estimator) float64 {
    return StdBlock(est)/math.Sqrt(float64(est.Nblock-1))
}

func PrintDetail(est Estimator) {
    mean := MeanBlock(est)
    err  := StdError(est)
    std  := Std(est)

    std2 := std*std
    tau  := err*err/std2*float64(est.N)
    tau   = (tau-1)*0.5

    fmt.Printf("------------------------ \n")
    fmt.Printf("Observable : %s \n",est.Name)
    fmt.Printf("Nblock=%v Bsize=%v \n",est.Nblock,est.Bsize)
    fmt.Printf("mean      = %e \n",mean)
    fmt.Printf("std error = %e \n",err)
    fmt.Printf("std div   = %e \n",std)
    fmt.Printf("cor time  = %e \n",tau)
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


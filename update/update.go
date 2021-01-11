package main

import (
    "fmt"
    "math"
    "math/rand"
)

func pointSamplingWeighted(t float64, w []float64) int {
    dis := rand.Float64()*t

    for i:=0;i<len(w);i++ {
        if dis<w[i] {
            return i
        }
        dis -= w[i]
    }
    fmt.Printf("pointSamplingWeighted : out of range!\n")
    return len(w)
}

var MaxKeyFromSample uint64 = 1000000000000000 // 10E15
func squenceSamplingUniform(lambda float64) []uint64{
    var sequence []uint64
    var k float64 = 0

    dis := rand.Float64()
    k -= math.Log(dis)/lambda
    for k<1.0 {
        sequence = append(sequence,uint64(k*float64(MaxKeyFromSample)))
        dis = rand.Float64()
        k -= math.Log(dis)/lambda
    }

    return sequence
}

func main() {
    rand.Seed(89374)
    w := []float64{1,2,3,4,5,6,7,8,9}
    n := []float64{0,0,0,0,0,0,0,0,0}

    t := float64(0)
    for i:=0;i<len(w);i++ {
        t += w[i]
    }

    nn := 1000
    for i:=0;i<nn;i++ {
        d := pointSamplingWeighted(t,w)
        n[d] += 1.0/float64(nn)
    }

    fmt.Println(n)
    fmt.Println(squenceSamplingUniform(10000.0))
}

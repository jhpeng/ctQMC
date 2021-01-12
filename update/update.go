package main

import (
    "fmt"
    "math"
    "math/rand"
    . "github.com/jhpeng/ctQMC/dtype"
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

func Remove(w WorldLine) {
    var key uint64
    var v   Vertex

    sequence1 := w.SequenceB
    sequence2 := w.SequenceA
    if w.Flag {
        sequence1 = w.SequenceA
        sequence2 = w.SequenceB
    }

    k:=0
    for i:=0;i<w.Nvertices;i++ {
        key = sequence1[i]
        v   = w.Table[key]
        check_delete := true

        for j:=0;j<v.HNspin;j++ {
            if v.State[j]!=v.State[j+v.HNspin]{
                check_delete = false
            }
        }

        if check_delete {
            delete(w.Table,key)
        } else {
            sequence2[k] = sequence1[i]
            k++
        }
    }
    w.Nvertices = k
    w.Flag = !(w.Flag)
}

func Insert(w WorldLine, m Model) {
    lambda := (m.Tweight)*(w.Beta)
    insert_seq := squenceSamplingUniform(lambda)

    length := len(insert_seq)+w.Nvertices+1024
    if length>len(w.SequenceA) {
        temp1 := make([]uint64,length)
        copy(temp1,w.SequenceA)
        temp2 := make([]uint64,length)
        copy(temp2,w.SequenceB)
        w.SequenceA = temp1
        w.SequenceB = temp2
    }

    sequence1 := w.SequenceB
    sequence2 := w.SequenceA
    if w.Flag {
        sequence1 = w.SequenceA
        sequence2 = w.SequenceB
    }

    pstate := make([]int,m.Nsite)
    copy(pstate,w.State)

    k:=0
    n:=0
    var v Vertex
    for i:=0;i<len(insert_seq);i++ {
        key1 := sequence1[k]
        key2 := insert_seq[i]

        for (key1<key2) && (k<w.Nvertices) {
            v = w.Table[key1]
            indices := m.Bond2index[v.Bond]
            for i_site:=0;i_site<v.HNspin;i_site++ {
                index := indices[i_site]
                pstate[index] = v.State[v.HNspin+i_site]
            }

            sequence2[n] = key1
            n++
            k++
            key1 = sequence1[k]
        }

        if key1!=key2 {
            bond := pointSamplingWeighted(m.Tweight,m.Bond2weight)
            t := m.Bond2type[bond]
            indices := m.Bond2index[bond]
            hNspin := m.Bond2hNspin[bond]
            rule := m.InsertRule[t]
            lstate := make([]int,2*hNspin)
            for i_site:=0;i_site<hNspin;i_site++ {
                lstate[i_site]        = pstate[indices[i_site]]
                lstate[i_site+hNspin] = pstate[indices[i_site]]
            }

            if rule(lstate) {
                var iv Vertex
                iv.Bond = bond
                iv.HNspin = hNspin
                iv.State = lstate

                w.Table[key2] = iv
                sequence2[n] = key2
                n++
            }
        }
    }

    for k<w.Nvertices {
        key1 := sequence1[k]
        sequence2[n] = key1
        n++
        k++
    }

    w.Nvertices = n
    w.Flag = !(w.Flag)
}

func main() {
    a := make([]int,65)
    fmt.Println(len(a),cap(a))
}

package update

import (
    "fmt"
    "math"
    "math/rand"
    "github.com/jhpeng/ctQMC/update/cluster"
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

var cmf []float64
var count []uint64
func pointSamplingWeighted2(t float64, w []float64) int {
    dis := rand.Float64()*t

    if cmf==nil {
        cmf    = make([]float64, len(w))
        count  = make([]uint64,  len(w))
        cmf[0] = w[0]
        count[0] = 0
        for i:=1;i<len(w);i++ {
            cmf[i] = cmf[i-1]+w[i]
            count[i] = 0
        }
    }

    i:=-1
    j:=len(cmf)
    d:=j-i
    for d>1 {
        if dis > cmf[i+d/2] {
            i = i+d/2
        } else {
            j = i+d/2
        }
        d = j-i
    }

    count[j]++

    return j
}

var MaxKeyFromSample uint64 = 1000000000000000 // 10E15
func squenceSamplingUniform(lambda float64, start uint64) []uint64{
    var sequence []uint64
    var k float64 = 0

    dis := rand.Float64()
    k -= math.Log(dis)/lambda
    for k<1.0 {
        sequence = append(sequence,uint64(k*float64(MaxKeyFromSample))+start)
        dis = rand.Float64()
        k -= math.Log(dis)/lambda
    }

    return sequence
}

func Remove(w WorldLine) WorldLine{
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

    return w
}

func Insert(w WorldLine, m Model) WorldLine{
    lambda := (m.Tweight)*(w.Beta)
    insert_seq := squenceSamplingUniform(lambda,0)

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
            bond := pointSamplingWeighted2(m.Tweight,m.Bond2weight)
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

    return w
}

func InnerLink(w WorldLine, m Model) (map[Id]Id, map[Id]int) {
    p := make(map[Id]Id)
    c := make(map[Id]int)

    sequence := w.SequenceB
    if w.Flag {
        sequence = w.SequenceA
    }

    for i:=0;i<w.Nvertices;i++ {
        key    := sequence[i]
        v      := w.Table[key]
        bond   := v.Bond
        hNspin := v.HNspin
        t      := m.Bond2type[bond]
        rule   := m.LinkRule[t]

        for j:=0;j<2*hNspin;j++ {
            idn := Id{Key:key, I:j}
            idv := Id{Key:key, I:rule[j]}
            p[idn] = idv
            c[idn] = rule[j+2*hNspin]
        }
    }

    return p,c
}

var null Id = Id{Key:0, I:-1}

func OuterLink(w WorldLine, m Model, p map[Id]Id, c map[Id]int) {
    last  := w.Last
    first := w.First

    for i:=0;i<w.Nsite;i++ {
        last[i] = null
        first[i] = null
    }

    sequence := w.SequenceB
    if w.Flag {
        sequence = w.SequenceA
    }

    for i:=0;i<w.Nvertices;i++ {
        key     := sequence[i]
        v       := w.Table[key]
        bond    := v.Bond
        hNspin  := v.HNspin
        indices := m.Bond2index[bond]

        for j:=0;j<hNspin;j++ {
            index := indices[j]
            idp := Id{Key:key, I:j}
            idn := Id{Key:key, I:j+hNspin}
            if first[index]==null {
                first[index] = idp
                last[index]  = idn
            } else {
                cluster.Union(p,c,last[index],idp)
                last[index] = idn
            }
        }
    }

    for i:=0;i<w.Nsite;i++ {
        if first[i]!=null {
            cluster.Union(p,c,first[i],last[i])
        }
    }
}

func FlipCluster(w WorldLine, p map[Id]Id) {
    newState := make(map[Id]int)
    first    := w.First

    sequence := w.SequenceB
    if w.Flag {
        sequence = w.SequenceA
    }

    for i:=0;i<w.Nvertices;i++ {
        key     := sequence[i]
        v       := w.Table[key]
        hNspin  := v.HNspin
        state   := v.State

        for j:=0;j<2*hNspin;j++ {
            idv := Id{Key:key, I:j}
            idr := cluster.Root(p,idv)
            if newState[idr]==0 {
                newState[idr] = 1
                if rand.Float64()<0.5 {
                    newState[idr]=-1
                }
            }

            state[j] = state[j]*newState[idr]
        }
    }

    for i:=0;i<w.Nsite;i++ {
        id  := first[i]
        if id!=null {
            key := id.Key
            j   := id.I

            v := w.Table[key]
            w.State[i] = v.State[j]
        } else if rand.Float64()<0.5 {
            w.State[i] =  1
        } else {
            w.State[i] = -1
        }
    }
}

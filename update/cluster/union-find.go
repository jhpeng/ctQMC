package cluster

func Root(p []int, v int) int {
    u := v
    for p[u]!=u {
        u = p[u]
    }

    return u
}

func Compress(p []int, v int, r int) {
    u := v
    for p[u]!=u {
        s := p[u]
        p[u] = r
        u = s
    }
}

func Union(p []int, w []int, va int, vb int){
    ra := Root(p,va)
    rb := Root(p,vb)
    r0 := ra
    r1 := rb
    if w[ra]<w[rb] {
        r0 = rb
        r1 = ra
    }

    w[r0] = w[r0]+w[r1]
    p[r1] = r0

    Compress(p,va,r0)
    Compress(p,vb,r0)
}


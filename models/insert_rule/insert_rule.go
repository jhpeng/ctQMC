package insert_rule


func InsertRuleIsing1(state []int) bool {
    return true
}

func InsertRuleIsing2(state []int) bool {
    if state[0]*state[1]==1{
        return true
    } else {
        return false
    }
}

func InsertRuleHeisenbergAFM(state []int) bool {
    if state[0]*state[1]==-1{
        return true
    } else {
        return false
    }
}

func InsertRuleJQ3(state []int) bool {
    if state[0]*state[1]==-1 {
        if state[2]*state[3]==-1 {
            if state[4]*state[5]==-1 {
                return true
            } else {
                return false
            }
        } else {
            return false
        }
    } else {
        return false
    }
}


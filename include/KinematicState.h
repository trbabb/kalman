/* 
 * File:   KinematicState.h
 * Author: tbabb
 *
 * Created on February 3, 2015, 3:13 PM
 */

#ifndef KINEMATICSTATE_H
#define	KINEMATICSTATE_H


// position
#define KINSTATE_IDX_X       0
// velocity
#define KINSTATE_IDX_V       3
// acceleration
#define KINSTATE_IDX_A       6
// angular velocity
#define KINSTATE_IDX_OMEGA   9
// absolute orientation
#define KINSTATE_IDX_ORIENT 12
// n elements
#define KINSTATE_SIZE       16


/*******************************************
 * State class                             *
 *******************************************/


template <typename T>
struct KinematicState {
    
    KinematicState():orient(0,0,0,1) {}
    
    Vec<T,3> x;
    Vec<T,3> v;
    Vec<T,3> a;
    Vec<T,3> omega;
    Quat<T>  orient;
    
    inline       T* begin()       { return x.begin();     }
    inline const T* begin() const { return x.begin();     }
    inline       T* end()         { return orient.end();  }
    inline const T* end()   const { return orient.end();  }
    inline index_t size()   const { return KINSTATE_SIZE; }
    
    KinematicState<T>& operator+=(const KinematicState &k) {
        x += k.x;
        v += k.v;
        a += k.a;
        omega  += k.omega;
        orient += k.orient;
    }
    
};


#endif	/* KINEMATICSTATE_H */


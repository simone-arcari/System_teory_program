function J = costIndex(w)
% COSTINDEX = This fuction returns the value J which is the cost index of control u

    w_t = transpose(w);
    J = w_t * w;





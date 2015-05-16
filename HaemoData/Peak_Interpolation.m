function LVP = Peak_Interpolation(LVP, pre_discont, post_discont)
% This function bridges the discontinuity at the peak of LVP trace by
% interpolating a quadratic function between two anchor points. 
% 
%% Located anchor points
if (LVP(pre_discont) > LVP(post_discont))
    post_anchor = post_discont + 1;
    pre_anchor = pre_discont;
    inter_pnt = pre_discont + 1;
else 
    pre_anchor = pre_discont - 1;
    post_anchor = post_discont;
    inter_pnt = post_discont - 1;
end
% 
LVP(inter_pnt) = (LVP(pre_anchor) + LVP(post_anchor))/2;
% %% Estimate gradients at each anchor point. 
% pre_grad = LVP(pre_anchor) - LVP(pre_anchor - 1);
% post_grad = LVP(post_anchor + 1) - LVP(post_anchor);
% 
% %% Calculate coefficients of cubic function
% a = 2*LVP(pre_anchor) - 2*LVP(post_anchor) + pre_grad + post_grad;
% b = 3*LVP(post_anchor) - 3*LVP(pre_anchor) - post_grad - 2*pre_grad;
% c = pre_grad;
% d = LVP(pre_anchor);
% 
% %% Calculate LVP of interpolated point
% LVP(inter_pnt) = a*(0.5^3) + b*0.5^2 + c*0.5 + d;

return
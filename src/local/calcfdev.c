#include <math.h>
#include "typedefs.h"
#include "main.h"

static real dF2dC12H(real c12H,real c12S,real c6S,real kbO,real kbH,rvec r[],
		     real qO,real qH) 
{
  /* PvM, Wed Nov 25 18:12:03 CET 1998
   *
   * This file contains output generated with the following maple command:
   *
   * C(diff(ftot2,c12H),optimized,filename=`yaw2_c12H.c`);
   *
   * for details see: yaw2_analytic.map 
   * 
   * This is the isotropic model!
   */
  real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
  real t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
  real t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
  real t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
  real t41, t42, t43, t44, t45, t46, t47, t48, t49, t50;
  real t51, t52, t53, t54, t55, t56, t57, t58, t59, t60;
  real t61, t62, t63, t64, t65, t66, t67, t68, t69, t70;
  real t71, t72, t73, t74, t75, t76, t77, t78, t79, t80;
  real t81, t82, t83, t84, t85, t86, t87, t88, t89, t90;
  real t91, t92, t93, t94, t95, t96, t97, t98, t99, t100;
  real t101, t102, t103, t104, t105, t106, t107, t108, t109, t110;
  real t111, t112, t113, t114, t115, t116, t117, t118, t119, t120;
  real t121, t122, t123, t124, t125, t126, t127, t128, t129, t130;
  real t131, t132, t133, t134, t135, t136, t137, t138, t139, t140;
  real t141, t142, t143, t144, t145, t146, t147, t148, t149, t150;
  real t151, t152, t153, t154, t155, t156, t157, t158, t159, t160;
  real t161, t162, t163, t164, t165, t166, t167, t168, t169, t170;
  real t171, t172, t173, t174, t175, t176, t177, t178, t179, t180;
  real t181, t182, t183, t184, t185, t186, t187, t188, t189, t190;
  real t191, t192, t193, t194, t195, t196, t197, t198, t199, t200;
  real t201, t202, t203, t204, t205, t206, t207, t208, t209, t210;
  real t211, t212, t213, t214, t215, t216, t217, t218, t219, t220;
  real t221, t222, t223, t224, t225, t226, t227, t228, t229, t230;
  real t231, t232, t233, t234, t235, t236, t237, t238, t239, t240;
  real t241, t242, t243, t244, t245, t246, t247, t248, t249, t250;
  real t251, t252, t253, t254, t255, t256, t257, t258, t259, t260;
  real t261, t262, t263, t264, t265, t266, t267, t268, t269, t270;
  real t271, t272, t273, t274, t275, t276, t277, t278, t279, t280;
  real t281, t282, t283, t284, t285, t286, t287, t288, t289, t290;
  real t291, t292, t293, t294, t295, t296, t297, t298, t299, t300;
  real t301, t302, t303, t304, t305, t306, t307, t308, t309, t310;
  real t311, t312, t313, t314, t315, t316, t317, t318, t319, t320;
  real t321, t322, t323, t324, t325, t326, t327, t328, t329, t330;
  real t331, t332, t333, t334, t335, t336, t337, t338, t339, t340;
  real t341, t342, t343, t344, t345, t346, t347, t348, t349, t350;
  real t351, t352, t353, t354, t355, t356, t357, t358, t359, t360;
  real t361, t362, t363;
  
  t1 = c12H*c12H;
  t2 = sqrt(t1);
  t4 = pow(r[2][XX]-r[6][XX],2.0);
  t6 = pow(r[2][YY]-r[6][YY],2.0);
  t8 = pow(r[2][ZZ]-r[6][ZZ],2.0);
  t9 = t4+t6+t8;
  t10 = t9*t9;
  t11 = t10*t10;
  t13 = 1/t11/t10;
  t15 = pow(r[2][XX],2.0);
  t17 = pow(r[6][XX],2.0);
  t18 = pow(r[2][YY],2.0);
  t20 = pow(r[6][YY],2.0);
  t21 = pow(r[2][ZZ],2.0);
  t23 = pow(r[6][ZZ],2.0);
  t25 = sqrt(t15-2.0*r[2][XX]*r[6][XX]+t17+t18-2.0*r[2][YY]*r[6][YY]+t20+t21
	     -2.0*r[2][ZZ]*r[6][ZZ]+t23);
  t26 = 1/t25;
  t29 = sqrt(c12H*c12S);
  t31 = pow(r[4][XX]-r[6][XX],2.0);
  t33 = pow(r[4][YY]-r[6][YY],2.0);
  t35 = pow(r[4][ZZ]-r[6][ZZ],2.0);
  t36 = t31+t33+t35;
  t37 = t36*t36;
  t38 = t37*t37;
  t40 = 1/t38/t37;
  t42 = pow(r[4][XX],2.0);
  t44 = pow(r[4][YY],2.0);
  t46 = pow(r[4][ZZ],2.0);
  t49 = sqrt(t42-2.0*r[4][XX]*r[6][XX]+t17+t44-2.0*r[4][YY]*r[6][YY]+t20+t46
	     -2.0*r[4][ZZ]*r[6][ZZ]+t23);
  t50 = 1/t49;
  t52 = qH*qH;
  t54 = pow(r[2][XX]-r[7][XX],2.0);
  t56 = pow(r[2][YY]-r[7][YY],2.0);
  t58 = pow(r[2][ZZ]-r[7][ZZ],2.0);
  t59 = t54+t56+t58;
  t62 = qO*qH;
  t64 = pow(r[1][XX]-r[6][XX],2.0);
  t66 = pow(r[1][YY]-r[6][YY],2.0);
  t68 = pow(r[1][ZZ]-r[6][ZZ],2.0);
  t72 = -qO-2.0*qH;
  t73 = qH*t72;
  t75 = pow(r[3][XX]-r[8][XX],2.0);
  t77 = pow(r[3][YY]-r[8][YY],2.0);
  t79 = pow(r[3][ZZ]-r[8][ZZ],2.0);
  t80 = t75+t77+t79;
  t84 = pow(r[3][XX]-r[5][XX],2.0);
  t86 = pow(r[3][YY]-r[5][YY],2.0);
  t88 = pow(r[3][ZZ]-r[5][ZZ],2.0);
  t95 = pow(r[2][XX]-r[8][XX],2.0);
  t97 = pow(r[2][YY]-r[8][YY],2.0);
  t99 = pow(r[2][ZZ]-r[8][ZZ],2.0);
  t100 = t95+t97+t99;
  t104 = pow(r[4][XX]-r[7][XX],2.0);
  t106 = pow(r[4][YY]-r[7][YY],2.0);
  t108 = pow(r[4][ZZ]-r[7][ZZ],2.0);
  t109 = t104+t106+t108;
  t113 = pow(r[2][XX]-r[5][XX],2.0);
  t115 = pow(r[2][YY]-r[5][YY],2.0);
  t117 = pow(r[2][ZZ]-r[5][ZZ],2.0);
  t121 = t80*t80;
  t122 = t121*t121;
  t124 = 1/t122/t121;
  t126 = pow(r[3][XX],2.0);
  t128 = pow(r[8][XX],2.0);
  t129 = pow(r[3][YY],2.0);
  t131 = pow(r[8][YY],2.0);
  t132 = pow(r[3][ZZ],2.0);
  t134 = pow(r[8][ZZ],2.0);
  t136 = sqrt(t126-2.0*r[3][XX]*r[8][XX]+t128+t129-2.0*r[3][YY]*r[8][YY]+t131+
	      t132-2.0*r[3][ZZ]*r[8][ZZ]+t134);
  t137 = 1/t136;
  t139 = t100*t100;
  t140 = t139*t139;
  t142 = 1/t140/t139;
  t148 = sqrt(t15-2.0*r[2][XX]*r[8][XX]+t128+t18-2.0*r[2][YY]*r[8][YY]+t131+t21
	      -2.0*r[2][ZZ]*r[8][ZZ]+t134);
  t149 = 1/t148;
  t151 = t72*qO;
  t153 = pow(r[4][XX]-r[5][XX],2.0);
  t155 = pow(r[4][YY]-r[5][YY],2.0);
  t157 = pow(r[4][ZZ]-r[5][ZZ],2.0);
  t162 = pow(r[3][XX]-r[7][XX],2.0);
  t164 = pow(r[3][YY]-r[7][YY],2.0);
  t166 = pow(r[3][ZZ]-r[7][ZZ],2.0);
  t167 = t162+t164+t166;
  t168 = t167*t167;
  t169 = t168*t168;
  t171 = 1/t169/t168;
  t174 = pow(r[7][XX],2.0);
  t176 = pow(r[7][YY],2.0);
  t178 = pow(r[7][ZZ],2.0);
  t180 = sqrt(t126-2.0*r[3][XX]*r[7][XX]+t174+t129-2.0*r[3][YY]*r[7][YY]+t176+
	      t132-2.0*r[3][ZZ]*r[7][ZZ]+t178);
  t181 = 1/t180;
  t185 = t59*t59;
  t186 = t185*t185;
  t188 = 1/t186/t185;
  t194 = sqrt(t15-2.0*r[2][XX]*r[7][XX]+t174+t18-2.0*r[2][YY]*r[7][YY]+t176+t21
	      -2.0*r[2][ZZ]*r[7][ZZ]+t178);
  t195 = 1/t194;
  t197 = 0.12E2*t2*t13*t26+0.12E2*t29*t40*t50+0.1389354102E3*t52/t59+
    0.1389354102E3*t62/(t64+t66+t68)+0.1389354102E3*t73/t80+0.1389354102E3*t62/(t84
										+t86+t88)+0.1389354102E3*t73/t36+0.1389354102E3*t73/t100+0.1389354102E3*t73/
    t109+0.1389354102E3*t62/(t113+t115+t117)+0.12E2*t29*t124*t137+0.12E2*t29*t142*
    t149+0.1389354102E3*t151/(t153+t155+t157)+0.12E2*t2*t171*t181+0.1389354102E3*
    t52/t9+0.12E2*t2*t188*t195;
  t198 = t109*t109;
  t199 = t198*t198;
  t201 = 1/t199/t198;
  t207 = sqrt(t42-2.0*r[4][XX]*r[7][XX]+t174+t44-2.0*r[4][YY]*r[7][YY]+t176+t46
	      -2.0*r[4][ZZ]*r[7][ZZ]+t178);
  t208 = 1/t207;
  t210 = c12S*c12S;
  t211 = sqrt(t210);
  t213 = pow(r[4][XX]-r[8][XX],2.0);
  t215 = pow(r[4][YY]-r[8][YY],2.0);
  t217 = pow(r[4][ZZ]-r[8][ZZ],2.0);
  t218 = t213+t215+t217;
  t219 = t218*t218;
  t220 = t219*t219;
  t228 = sqrt(t42-2.0*r[4][XX]*r[8][XX]+t128+t44-2.0*r[4][YY]*r[8][YY]+t131+t46
	      -2.0*r[4][ZZ]*r[8][ZZ]+t134);
  t229 = 1/t228;
  t231 = c6S*c6S;
  t232 = sqrt(t231);
  t238 = pow(r[3][XX]-r[6][XX],2.0);
  t240 = pow(r[3][YY]-r[6][YY],2.0);
  t242 = pow(r[3][ZZ]-r[6][ZZ],2.0);
  t243 = t238+t240+t242;
  t244 = t243*t243;
  t245 = t244*t244;
  t247 = 1/t245/t244;
  t253 = sqrt(t126-2.0*r[3][XX]*r[6][XX]+t17+t129-2.0*r[3][YY]*r[6][YY]+t20+
	      t132-2.0*r[3][ZZ]*r[6][ZZ]+t23);
  t254 = 1/t253;
  t257 = pow(r[1][XX]-r[7][XX],2.0);
  t259 = pow(r[1][YY]-r[7][YY],2.0);
  t261 = pow(r[1][ZZ]-r[7][ZZ],2.0);
  t266 = pow(r[1][XX]-r[8][XX],2.0);
  t268 = pow(r[1][YY]-r[8][YY],2.0);
  t270 = pow(r[1][ZZ]-r[8][ZZ],2.0);
  t277 = pow(r[5][XX],2.0);
  t279 = pow(r[5][YY],2.0);
  t281 = pow(r[5][ZZ],2.0);
  t283 = sqrt(t128-2.0*r[8][XX]*r[5][XX]+t277+t131-2.0*r[8][YY]*r[5][YY]+t279+
	      t134-2.0*r[8][ZZ]*r[5][ZZ]+t281);
  t285 = t72*t72;
  t294 = sqrt(t128-2.0*r[8][XX]*r[6][XX]+t17+t131-2.0*r[8][YY]*r[6][YY]+t20+
	      t134-2.0*r[8][ZZ]*r[6][ZZ]+t23);
  t300 = sqrt(t42-2.0*r[4][XX]*r[2][XX]+t15+t44-2.0*r[4][YY]*r[2][YY]+t18+t46
	      -2.0*r[4][ZZ]*r[2][ZZ]+t21);
  t306 = sqrt(t42-2.0*r[4][XX]*r[3][XX]+t126+t44-2.0*r[4][YY]*r[3][YY]+t129+t46
	      -2.0*r[4][ZZ]*r[3][ZZ]+t132);
  t308 = qO*qO;
  t310 = pow(r[1][XX]-r[5][XX],2.0);
  t312 = pow(r[1][YY]-r[5][YY],2.0);
  t314 = pow(r[1][ZZ]-r[5][ZZ],2.0);
  t319 = pow(r[1][XX],2.0);
  t321 = pow(r[1][YY],2.0);
  t323 = pow(r[1][ZZ],2.0);
  t325 = sqrt(t42-2.0*r[4][XX]*r[1][XX]+t319+t44-2.0*r[4][YY]*r[1][YY]+t321+t46
	      -2.0*r[4][ZZ]*r[1][ZZ]+t323);
  t331 = sqrt(t128-2.0*r[8][XX]*r[7][XX]+t174+t131-2.0*r[8][YY]*r[7][YY]+t176+
	      t134-2.0*r[8][ZZ]*r[7][ZZ]+t178);
  t333 = 0.12E2*t29*t201*t208+0.12E2*t211/t220/t219*t229-0.6E1*t232/t219/
    t218*t229+0.12E2*t2*t247*t254+0.1389354102E3*t62/(t257+t259+t261)+
    0.1389354102E3*t151/(t266+t268+t270)+0.1389354102E3*t52/t167-kbO*t283+
    0.1389354102E3*t285/t218+0.1389354102E3*t52/t243-kbH*t294-kbH*t300-kbH*t306+
    0.1389354102E3*t308/(t310+t312+t314)-kbO*t325-kbH*t331;
  t335 = 1/t29;
  t339 = 1/t2;
  t363 = 2.0*(t197+t333)*(0.6E1*t335*t201*t208*c12S+0.12E2*t339*t171*t181*
			  c12H+0.6E1*t335*t40*t50*c12S+0.6E1*t335*t124*t137*c12S+0.12E2*t339*t247*t254*
			  c12H+0.6E1*t335*t142*t149*c12S+0.12E2*t339*t188*t195*c12H+0.12E2*t339*t13*t26*
			  c12H);
			  
  return t363;
}

static real dF2dC12S(real c12H,real c12S,real c6S,real kbO,real kbH,rvec r[],
		     real qO,real qH) 
{
  /* PvM, Wed Nov 25 18:12:03 CET 1998
   *
   * This file contains output generated with the following maple command:
   *
   * C(diff(ftot2,c12S),optimized,filename=`yaw2_c12S.c`);
   *
   * for details see: yaw2_analytic.map 
   * 
   */
   
  real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
  real t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
  real t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
  real t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
  real t41, t42, t43, t44, t45, t46, t47, t48, t49, t50;
  real t51, t52, t53, t54, t55, t56, t57, t58, t59, t60;
  real t61, t62, t63, t64, t65, t66, t67, t68, t69, t70;
  real t71, t72, t73, t74, t75, t76, t77, t78, t79, t80;
  real t81, t82, t83, t84, t85, t86, t87, t88, t89, t90;
  real t91, t92, t93, t94, t95, t96, t97, t98, t99, t100;
  real t101, t102, t103, t104, t105, t106, t107, t108, t109, t110;
  real t111, t112, t113, t114, t115, t116, t117, t118, t119, t120;
  real t121, t122, t123, t124, t125, t126, t127, t128, t129, t130;
  real t131, t132, t133, t134, t135, t136, t137, t138, t139, t140;
  real t141, t142, t143, t144, t145, t146, t147, t148, t149, t150;
  real t151, t152, t153, t154, t155, t156, t157, t158, t159, t160;
  real t161, t162, t163, t164, t165, t166, t167, t168, t169, t170;
  real t171, t172, t173, t174, t175, t176, t177, t178, t179, t180;
  real t181, t182, t183, t184, t185, t186, t187, t188, t189, t190;
  real t191, t192, t193, t194, t195, t196, t197, t198, t199, t200;
  real t201, t202, t203, t204, t205, t206, t207, t208, t209, t210;
  real t211, t212, t213, t214, t215, t216, t217, t218, t219, t220;
  real t221, t222, t223, t224, t225, t226, t227, t228, t229, t230;
  real t231, t232, t233, t234, t235, t236, t237, t238, t239, t240;
  real t241, t242, t243, t244, t245, t246, t247, t248, t249, t250;
  real t251, t252, t253, t254, t255, t256, t257, t258, t259, t260;
  real t261, t262, t263, t264, t265, t266, t267, t268, t269, t270;
  real t271, t272, t273, t274, t275, t276, t277, t278, t279, t280;
  real t281, t282, t283, t284, t285, t286, t287, t288, t289, t290;
  real t291, t292, t293, t294, t295, t296, t297, t298, t299, t300;
  real t301, t302, t303, t304, t305, t306, t307, t308, t309, t310;
  real t311, t312, t313, t314, t315, t316, t317, t318, t319, t320;
  real t321, t322, t323, t324, t325, t326, t327, t328, t329, t330;
  real t331, t332, t333, t334, t335, t336, t337, t338, t339, t340;
  real t341, t342, t343, t344, t345, t346, t347, t348, t349, t350;
  real t351, t352, t353, t354, t355, t356, t357, t358, t359, t360;
   
  t1 = -qO-2.0*qH;
  t2 = qH*t1;
  t4 = pow(r[3][XX]-r[8][XX],2.0);
  t6 = pow(r[3][YY]-r[8][YY],2.0);
  t8 = pow(r[3][ZZ]-r[8][ZZ],2.0);
  t9 = t4+t6+t8;
  t13 = sqrt(c12H*c12S);
  t15 = pow(r[4][XX]-r[6][XX],2.0);
  t17 = pow(r[4][YY]-r[6][YY],2.0);
  t19 = pow(r[4][ZZ]-r[6][ZZ],2.0);
  t20 = t15+t17+t19;
  t21 = t20*t20;
  t22 = t21*t21;
  t24 = 1/t22/t21;
  t26 = pow(r[4][XX],2.0);
  t28 = pow(r[6][XX],2.0);
  t29 = pow(r[4][YY],2.0);
  t31 = pow(r[6][YY],2.0);
  t32 = pow(r[4][ZZ],2.0);
  t34 = pow(r[6][ZZ],2.0);
  t36 = sqrt(t26-2.0*r[4][XX]*r[6][XX]+t28+t29-2.0*r[4][YY]*r[6][YY]+t31+t32
	     -2.0*r[4][ZZ]*r[6][ZZ]+t34);
  t37 = 1/t36;
  t39 = c12H*c12H;
  t40 = sqrt(t39);
  t42 = pow(r[3][XX]-r[7][XX],2.0);
  t44 = pow(r[3][YY]-r[7][YY],2.0);
  t46 = pow(r[3][ZZ]-r[7][ZZ],2.0);
  t47 = t42+t44+t46;
  t48 = t47*t47;
  t49 = t48*t48;
  t53 = pow(r[3][XX],2.0);
  t55 = pow(r[7][XX],2.0);
  t56 = pow(r[3][YY],2.0);
  t58 = pow(r[7][YY],2.0);
  t59 = pow(r[3][ZZ],2.0);
  t61 = pow(r[7][ZZ],2.0);
  t63 = sqrt(t53-2.0*r[3][XX]*r[7][XX]+t55+t56-2.0*r[3][YY]*r[7][YY]+t58+t59
	     -2.0*r[3][ZZ]*r[7][ZZ]+t61);
  t67 = pow(r[3][XX]-r[6][XX],2.0);
  t69 = pow(r[3][YY]-r[6][YY],2.0);
  t71 = pow(r[3][ZZ]-r[6][ZZ],2.0);
  t72 = t67+t69+t71;
  t73 = t72*t72;
  t74 = t73*t73;
  t82 = sqrt(t53-2.0*r[3][XX]*r[6][XX]+t28+t56-2.0*r[3][YY]*r[6][YY]+t31+t59
	     -2.0*r[3][ZZ]*r[6][ZZ]+t34);
  t86 = pow(r[2][XX]-r[7][XX],2.0);
  t88 = pow(r[2][YY]-r[7][YY],2.0);
  t90 = pow(r[2][ZZ]-r[7][ZZ],2.0);
  t91 = t86+t88+t90;
  t92 = t91*t91;
  t93 = t92*t92;
  t97 = pow(r[2][XX],2.0);
  t99 = pow(r[2][YY],2.0);
  t101 = pow(r[2][ZZ],2.0);
  t104 = sqrt(t97-2.0*r[2][XX]*r[7][XX]+t55+t99-2.0*r[2][YY]*r[7][YY]+t58+t101
	      -2.0*r[2][ZZ]*r[7][ZZ]+t61);
  t108 = pow(r[4][XX]-r[7][XX],2.0);
  t110 = pow(r[4][YY]-r[7][YY],2.0);
  t112 = pow(r[4][ZZ]-r[7][ZZ],2.0);
  t113 = t108+t110+t112;
  t116 = t113*t113;
  t117 = t116*t116;
  t119 = 1/t117/t116;
  t125 = sqrt(t26-2.0*r[4][XX]*r[7][XX]+t55+t29-2.0*r[4][YY]*r[7][YY]+t58+t32
	      -2.0*r[4][ZZ]*r[7][ZZ]+t61);
  t126 = 1/t125;
  t128 = qH*qO;
  t130 = pow(r[2][XX]-r[5][XX],2.0);
  t132 = pow(r[2][YY]-r[5][YY],2.0);
  t134 = pow(r[2][ZZ]-r[5][ZZ],2.0);
  t138 = c6S*c6S;
  t139 = sqrt(t138);
  t141 = pow(r[4][XX]-r[8][XX],2.0);
  t143 = pow(r[4][YY]-r[8][YY],2.0);
  t145 = pow(r[4][ZZ]-r[8][ZZ],2.0);
  t146 = t141+t143+t145;
  t147 = t146*t146;
  t152 = pow(r[8][XX],2.0);
  t154 = pow(r[8][YY],2.0);
  t156 = pow(r[8][ZZ],2.0);
  t158 = sqrt(t26-2.0*r[4][XX]*r[8][XX]+t152+t29-2.0*r[4][YY]*r[8][YY]+t154+t32
	      -2.0*r[4][ZZ]*r[8][ZZ]+t156);
  t159 = 1/t158;
  t162 = pow(r[1][XX]-r[6][XX],2.0);
  t164 = pow(r[1][YY]-r[6][YY],2.0);
  t166 = pow(r[1][ZZ]-r[6][ZZ],2.0);
  t170 = c12S*c12S;
  t171 = sqrt(t170);
  t172 = t147*t147;
  t174 = 1/t172/t147;
  t177 = t1*qO;
  t179 = pow(r[4][XX]-r[5][XX],2.0);
  t181 = pow(r[4][YY]-r[5][YY],2.0);
  t183 = pow(r[4][ZZ]-r[5][ZZ],2.0);
  t188 = pow(r[1][XX]-r[7][XX],2.0);
  t190 = pow(r[1][YY]-r[7][YY],2.0);
  t192 = pow(r[1][ZZ]-r[7][ZZ],2.0);
  t197 = pow(r[1][XX]-r[8][XX],2.0);
  t199 = pow(r[1][YY]-r[8][YY],2.0);
  t201 = pow(r[1][ZZ]-r[8][ZZ],2.0);
  t207 = t9*t9;
  t208 = t207*t207;
  t210 = 1/t208/t207;
  t216 = sqrt(t53-2.0*r[3][XX]*r[8][XX]+t152+t56-2.0*r[3][YY]*r[8][YY]+t154+t59
	      -2.0*r[3][ZZ]*r[8][ZZ]+t156);
  t217 = 1/t216;
  t219 = 0.1389354102E3*t2/t9+0.12E2*t13*t24*t37+0.12E2*t40/t49/t48/t63+
    0.12E2*t40/t74/t73/t82+0.12E2*t40/t93/t92/t104+0.1389354102E3*t2/t113+0.12E2*
    t13*t119*t126+0.1389354102E3*t128/(t130+t132+t134)-0.6E1*t139/t147/t146*t159+
    0.1389354102E3*t128/(t162+t164+t166)+0.12E2*t171*t174*t159+0.1389354102E3*t177/
    (t179+t181+t183)+0.1389354102E3*t128/(t188+t190+t192)+0.1389354102E3*t177/(t197
									       +t199+t201)+0.1389354102E3*t2/t20+0.12E2*t13*t210*t217;
  t221 = pow(r[3][XX]-r[5][XX],2.0);
  t223 = pow(r[3][YY]-r[5][YY],2.0);
  t225 = pow(r[3][ZZ]-r[5][ZZ],2.0);
  t230 = pow(r[2][XX]-r[6][XX],2.0);
  t232 = pow(r[2][YY]-r[6][YY],2.0);
  t234 = pow(r[2][ZZ]-r[6][ZZ],2.0);
  t235 = t230+t232+t234;
  t236 = t235*t235;
  t237 = t236*t236;
  t245 = sqrt(t97-2.0*r[2][XX]*r[6][XX]+t28+t99-2.0*r[2][YY]*r[6][YY]+t31+t101
	      -2.0*r[2][ZZ]*r[6][ZZ]+t34);
  t249 = pow(r[2][XX]-r[8][XX],2.0);
  t251 = pow(r[2][YY]-r[8][YY],2.0);
  t253 = pow(r[2][ZZ]-r[8][ZZ],2.0);
  t254 = t249+t251+t253;
  t255 = t254*t254;
  t256 = t255*t255;
  t258 = 1/t256/t255;
  t264 = sqrt(t97-2.0*r[2][XX]*r[8][XX]+t152+t99-2.0*r[2][YY]*r[8][YY]+t154+
	      t101-2.0*r[2][ZZ]*r[8][ZZ]+t156);
  t265 = 1/t264;
  t269 = qH*qH;
  t276 = sqrt(t152-2.0*r[8][XX]*r[6][XX]+t28+t154-2.0*r[8][YY]*r[6][YY]+t31+
	      t156-2.0*r[8][ZZ]*r[6][ZZ]+t34);
  t280 = t1*t1;
  t284 = pow(r[5][XX],2.0);
  t286 = pow(r[5][YY],2.0);
  t288 = pow(r[5][ZZ],2.0);
  t290 = sqrt(t152-2.0*r[8][XX]*r[5][XX]+t284+t154-2.0*r[8][YY]*r[5][YY]+t286+
	      t156-2.0*r[8][ZZ]*r[5][ZZ]+t288);
  t296 = sqrt(t26-2.0*r[4][XX]*r[3][XX]+t53+t29-2.0*r[4][YY]*r[3][YY]+t56+t32
	      -2.0*r[4][ZZ]*r[3][ZZ]+t59);
  t302 = sqrt(t26-2.0*r[4][XX]*r[2][XX]+t97+t29-2.0*r[4][YY]*r[2][YY]+t99+t32
	      -2.0*r[4][ZZ]*r[2][ZZ]+t101);
  t305 = pow(r[1][XX],2.0);
  t307 = pow(r[1][YY],2.0);
  t309 = pow(r[1][ZZ],2.0);
  t311 = sqrt(t26-2.0*r[4][XX]*r[1][XX]+t305+t29-2.0*r[4][YY]*r[1][YY]+t307+t32
	      -2.0*r[4][ZZ]*r[1][ZZ]+t309);
  t313 = qO*qO;
  t315 = pow(r[1][XX]-r[5][XX],2.0);
  t317 = pow(r[1][YY]-r[5][YY],2.0);
  t319 = pow(r[1][ZZ]-r[5][ZZ],2.0);
  t331 = sqrt(t152-2.0*r[8][XX]*r[7][XX]+t55+t154-2.0*r[8][YY]*r[7][YY]+t58+
	      t156-2.0*r[8][ZZ]*r[7][ZZ]+t61);
  t333 = 0.1389354102E3*t128/(t221+t223+t225)+0.12E2*t40/t237/t236/t245+
    0.12E2*t13*t258*t265+0.1389354102E3*t2/t254+0.1389354102E3*t269/t47-kbH*t276+
    0.1389354102E3*t269/t72+0.1389354102E3*t280/t146-kbO*t290-kbH*t296-kbH*t302-kbO
    *t311+0.1389354102E3*t313/(t315+t317+t319)+0.1389354102E3*t269/t91+
    0.1389354102E3*t269/t235-kbH*t331;
  t335 = 1/t13;
  t354 = 2.0*(t219+t333)*(0.6E1*t335*t119*t126*c12H+0.6E1*t335*t24*t37*c12H
			  +0.6E1*t335*t210*t217*c12H+0.6E1*t335*t258*t265*c12H+0.12E2/t171*t174*t159*c12S
			  );
  
  return t354;
}

void calc_f_dev(int natoms,real charge[],rvec x[],t_idef *idef,real *xiH,real *xiS)
{
  real qO,qH,kbO,kbH,c6S,c12S,c12H;
  int  i;
  
  qO   = charge[0];
  qH   = charge[1];
  c12H = idef->iparams[4].lj.c12;
  c6S  = idef->iparams[8].lj.c6;
  c12S = idef->iparams[8].lj.c12;
  kbO  = idef->iparams[9].harmonic.krA;
  kbH  = idef->iparams[10].harmonic.krA;

  fprintf(stdlog,"qO=%g, qH=%g, c12H=%g, c6H=%g\n",qO,qH,c12H,c12S);
  
  *xiH = dF2dC12H(c12H,c12S,c6S,kbO,kbH,x-1,qO,qH);
  *xiS = dF2dC12S(c12H,c12S,c6S,kbO,kbH,x-1,qO,qH);
}

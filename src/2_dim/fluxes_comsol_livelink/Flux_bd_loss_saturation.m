function f = Flux(N1,N2,x1,y1,x2,y2,rE,alpha1,alpha2,gamma1,gamma2)

%---------------------------------------------------------------------------------------------

%This function computes the fluxes for the arrangements with N1 first and N2 second catalysts
%with the coordinates x1,y1 and x2,y2, the catalyst radius rE in units of the system,
%and alpha1,alpha2 the dimensionless reaction-diffusion parameters for the first and second catalyst.
%The loss mechanism for the intermediate is an absorbing boundary condition.
%the reaction kinetics is Michaelis Menten with saturation parameters gamma1 and gamma2.

%---------------------------------------------------------------------------------------------

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

for i =1:N1
    
model.param.set(sprintf('x1%i',i), x1(i), 'x1 coordinate');

model.param.set(sprintf('y1%i',i), y1(i), 'y1 coordinate');


end

for i =1:N2
    
model.param.set(sprintf('x2%i',i), x2(i), 'x2 coordinate');

model.param.set(sprintf('y2%i',i), y2(i), 'y2 coordinate');


end



model.param.set('rE', rE, 'Enzyme radius');
model.param.set('alpha1', alpha1, 'reaction rate 1');
model.param.set('alpha2', alpha2, 'reaction rate 2');
model.param.set('gamma1', gamma1, 'saturation parameter 1');
model.param.set('gamma2', gamma2, 'saturation parameter 2');
model.param.set('s0', 1, 'boundary substrate concentration');

model.modelNode.create('mod1');
model.modelNode('mod1').name('Model 1');

model.geom.create('geom1', 2);
model.geom('geom1').selection.create('csel1', 'CumulativeSelection');
model.geom('geom1').selection('csel1').name('E2');
model.geom('geom1').selection.create('csel2', 'CumulativeSelection');
model.geom('geom1').selection('csel2').name('E1');
model.geom('geom1').selection.create('csel3', 'CumulativeSelection');
model.geom('geom1').selection('csel3').name('system');

for k = 1:(N1+N2)+1

model.geom('geom1').feature.create(sprintf('c%i',k), 'Circle');

end

model.geom('geom1').feature.create('dif1', 'Difference');
model.geom('geom1').feature.create('dif2', 'Difference');



model.geom('geom1').feature(sprintf('c%i',(N1+N2)+1)).set('contributeto', 'csel3');
model.geom('geom1').feature(sprintf('c%i',(N1+N2)+1)).set('createselection', 'on');

for j = 1:N1
    
model.geom('geom1').feature(sprintf('c%i',j)).set('r', 'rE');    
model.geom('geom1').feature(sprintf('c%i',j)).set('pos', {sprintf('x1%i',j) sprintf('y1%i',j)});
model.geom('geom1').feature(sprintf('c%i',j)).set('contributeto', 'csel1');
model.geom('geom1').feature(sprintf('c%i',j)).set('createselection', 'on');

end

for l = N1+1:(N1+N2)

model.geom('geom1').feature(sprintf('c%i',l)).set('r', 'rE');    
model.geom('geom1').feature(sprintf('c%i',l)).set('pos', {sprintf('x2%i',l-N1) sprintf('y2%i',l-N1)});
model.geom('geom1').feature(sprintf('c%i',l)).set('contributeto', 'csel2');
model.geom('geom1').feature(sprintf('c%i',l)).set('createselection', 'on');

end



model.geom('geom1').feature('dif1').selection('input').named('csel3');
model.geom('geom1').feature('dif1').selection('input2').named('csel1');
model.geom('geom1').feature('dif2').selection('input').named('csel3');
model.geom('geom1').feature('dif2').selection('input2').named('csel2');
model.geom('geom1').run;

model.physics.create('c', 'CoefficientFormPDE', 'geom1');
model.physics('c').feature.create('flux1', 'FluxBoundary', 1);
model.physics('c').feature('flux1').selection.named('geom1_csel1_bnd');
model.physics('c').feature.create('flux2', 'FluxBoundary', 1);
model.physics('c').feature('flux2').selection.named('geom1_csel2_bnd');
model.physics('c').feature.create('dir1', 'DirichletBoundary', 1);
model.physics('c').feature('dir1').selection.named('geom1_csel3_bnd');
model.physics.create('c2', 'CoefficientFormPDE', 'geom1');
model.physics('c2').field('dimensionless').field('s');
model.physics('c2').field('dimensionless').component({'s'});
model.physics('c2').feature.create('flux1', 'FluxBoundary', 1);
model.physics('c2').feature('flux1').selection.named('geom1_csel1_bnd');
model.physics('c2').feature.create('dir1', 'DirichletBoundary', 1);
model.physics('c2').feature('dir1').selection.named('geom1_csel3_bnd');

model.mesh.create('mesh1', 'geom1');
model.mesh('mesh1').feature.create('ftri1', 'FreeTri');

model.physics('c').name('PDE');
model.physics('c').prop('ShapeProperty').set('boundaryFlux', '0');
model.physics('c').feature('cfeq1').set('da', '0');
model.physics('c').feature('cfeq1').set('f', '0');
model.physics('c').feature('flux1').set('g', 'alpha1*s/(s/gamma1+1)');
model.physics('c').feature('flux2').set('q', 'alpha2/(u/gamma2+1)');
model.physics('c2').name('PDE 2');
model.physics('c2').prop('ShapeProperty').set('boundaryFlux', '0');
model.physics('c2').feature('cfeq1').set('da', '0');
model.physics('c2').feature('cfeq1').set('f', '0');
model.physics('c2').feature('flux1').set('q', 'alpha1/(s/gamma1+1)');
model.physics('c2').feature('dir1').set('r', 's0');

model.mesh('mesh1').feature('size').set('hauto', 1);
model.mesh('mesh1').feature('size').set('custom', 'on');
model.mesh('mesh1').feature('size').set('hgrad', '1.5');
model.mesh('mesh1').feature('size').set('hcurve', '0.1');
model.mesh('mesh1').feature('size').set('hmax', '0.1');
model.mesh('mesh1').feature('size').set('hmin', '0.001');
model.mesh('mesh1').run;

model.study.create('std1');
model.study('std1').feature.create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature.create('s1', 'Stationary');
model.sol('sol1').feature('s1').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature.remove('fcDef');

model.study('std1').feature('stat').set('initstudyhide', 'on');
model.study('std1').feature('stat').set('initsolhide', 'on');
model.study('std1').feature('stat').set('notstudyhide', 'on');
model.study('std1').feature('stat').set('notsolhide', 'on');

model.result.numerical.create('int1', 'IntLine');
model.result.numerical.create('int2', 'IntLine');
model.result.numerical('int1').selection.named('geom1_csel1_bnd');
model.result.numerical('int1').set('probetag', 'none');
model.result.numerical('int2').selection.named('geom1_csel2_bnd');
model.result.numerical('int2').set('probetag', 'none');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').name('Compile Equations: Stationary');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').feature('s1').set('control', 'stat');
model.sol('sol1').feature('s1').feature('dDef').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('fc1').set('probesel', 'manual');
model.sol('sol1').feature('s1').feature('fc1').set('termonres', 'off');
model.sol('sol1').runAll;

model.result.numerical('int1').set('descr', 's*alpha1/(s/gamma1+1)');
model.result.numerical('int1').set('expr', 's*alpha1/(s/gamma1+1)');
model.result.numerical('int2').set('expr', 'u*alpha2/(u/gamma2+1)');
model.result.numerical('int2').set('descr', 'u*alpha2/(u/gamma2+1)');
model.result.numerical('int1').setResult;
model.result.numerical('int2').setResult;
f(1)=model.result.numerical('int1').getReal;
f(2)=model.result.numerical('int2').getReal;
end
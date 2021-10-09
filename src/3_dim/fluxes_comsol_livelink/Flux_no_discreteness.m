function J = Flux_no_discreteness(N1,N2,x1,y1,z1,x2,y2,z2,rE,alpha1,alpha2)
%---------------------------------------------------------------------------------------------

%This function computes the fluxes for the arrangements with N1 first and N2 second catalysts
%with the coordinates x1,y1 and x2,y2, the catalyst radius rE in units of the system. The catalysts
%are reactive throughout their volume.
%alpha1,alpha2 are the dimensionless reaction-diffusion parameters for the first and second catalyst.
%The loss mechanism for the intermediate is an absorbing boundary condition.

%---------------------------------------------------------------------------------------------

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

for i =1:N1
    
model.param.set(sprintf('x1%i',i), x1(i), 'x coordinate');

model.param.set(sprintf('y1%i',i), y1(i), 'y coordinate');

model.param.set(sprintf('z1%i',i), z1(i), 'z coordinate');

end

for i =1:N2
    
model.param.set(sprintf('x2%i',i), x2(i), 'x coordinate');

model.param.set(sprintf('y2%i',i), y2(i), 'y coordinate');

model.param.set(sprintf('z2%i',i), z2(i), 'z coordinate');

end

model.param.set('alpha1', alpha1, 'reaction rate 1');
model.param.set('alpha2', alpha2, 'reaction rate 2');
model.param.set('rE', rE, 'Enzyme radius');
model.param.set('s0', '1');

model.modelNode.create('comp1');
model.modelNode('comp1').defineLocalCoord(false);

model.geom.create('geom1', 3);
model.geom('geom1').geomRep('comsol');
model.geom('geom1').selection.create('csel1', 'CumulativeSelection');
model.geom('geom1').selection('csel1').name('system');
model.geom('geom1').selection.create('csel2', 'CumulativeSelection');
model.geom('geom1').selection('csel2').name('E1');
model.geom('geom1').selection.create('csel3', 'CumulativeSelection');
model.geom('geom1').selection('csel3').name('E2');

for k = 1:N1+N2+1

model.geom('geom1').feature.create(sprintf('sph%i',k), 'Sphere');

end

model.geom('geom1').feature(sprintf('sph%i',(N1+N2+1))).set('contributeto', 'csel1');
model.geom('geom1').feature(sprintf('sph%i',(N1+N2+1))).set('createselection', 'on');

for j = 1:N1

model.geom('geom1').feature(sprintf('sph%i',j)).set('r', 'rE');
model.geom('geom1').feature(sprintf('sph%i',j)).set('contributeto', 'csel2');
model.geom('geom1').feature(sprintf('sph%i',j)).set('createselection', 'on');
model.geom('geom1').feature(sprintf('sph%i',j)).set('pos', {sprintf('x1%i',j) sprintf('y1%i',j) sprintf('z1%i',j)});

end

for j = 1:N2

model.geom('geom1').feature(sprintf('sph%i',j+N1)).set('r', 'rE');
model.geom('geom1').feature(sprintf('sph%i',j+N1)).set('contributeto', 'csel3');
model.geom('geom1').feature(sprintf('sph%i',j+N1)).set('createselection', 'on');
model.geom('geom1').feature(sprintf('sph%i',j+N1)).set('pos', {sprintf('x2%i',j) sprintf('y2%i',j) sprintf('z2%i',j)});

end

model.geom('geom1').run;
model.geom('geom1').run('fin');

model.physics.create('c', 'CoefficientFormPDE', 'geom1');
model.physics('c').field('dimensionless').field('s');
model.physics('c').field('dimensionless').component({'s'});
model.physics('c').selection.named('geom1_csel1_dom');
model.physics('c').create('dir1', 'DirichletBoundary', 2);
model.physics('c').feature('dir1').selection.named('geom1_csel1_bnd');
model.physics('c').create('cfeq2', 'CoefficientFormPDE', 3);
model.physics('c').feature('cfeq2').selection.named('geom1_csel2_dom');
model.physics.create('c2', 'CoefficientFormPDE', 'geom1');
model.physics('c2').field('dimensionless').field('i');
model.physics('c2').field('dimensionless').component({'i'});
model.physics('c2').selection.named('geom1_csel1_dom');
model.physics('c2').create('dir1', 'DirichletBoundary', 2);
model.physics('c2').feature('dir1').selection.named('geom1_csel1_bnd');
model.physics('c2').create('cfeq2', 'CoefficientFormPDE', 3);
model.physics('c2').feature('cfeq2').selection.named('geom1_csel2_dom');
model.physics('c2').create('cfeq3', 'CoefficientFormPDE', 3);
model.physics('c2').feature('cfeq3').selection.named('geom1_csel3_dom');

model.mesh.create('mesh1', 'geom1');
model.mesh('mesh1').feature.create('ftet1', 'FreeTet');

model.physics('c').feature('cfeq1').set('f', '0');
model.physics('c').feature('cfeq1').set('da', '0');
model.physics('c').feature('dir1').set('r', 's0');
model.physics('c').feature('cfeq2').set('a', 'alpha1');
model.physics('c').feature('cfeq2').set('f', '0');
model.physics('c').feature('cfeq2').set('da', '0');
model.physics('c2').feature('cfeq1').set('f', '0');
model.physics('c2').feature('cfeq1').set('da', '0');
model.physics('c2').feature('cfeq2').set('f', 'alpha1*s');
model.physics('c2').feature('cfeq2').set('da', '0');
model.physics('c2').feature('cfeq3').set('a', 'alpha2');
model.physics('c2').feature('cfeq3').set('f', '0');
model.physics('c2').feature('cfeq3').set('da', '0');

model.mesh('mesh1').feature('size').set('hauto', 3);
model.mesh('mesh1').feature('size').set('custom', 'on');
model.mesh('mesh1').feature('size').set('hgrad', '2');
model.mesh('mesh1').feature('size').set('hcurve', '0.1');
model.mesh('mesh1').feature('size').set('hmin', '0.001');
model.mesh('mesh1').feature('size').set('hmax', '0.5');
model.mesh('mesh1').run;

model.study.create('std1');
model.study('std1').feature.create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature.create('s1', 'Stationary');
model.sol('sol1').feature('s1').feature.create('se1', 'Segregated');
model.sol('sol1').feature('s1').feature('se1').feature.create('ss1', 'SegregatedStep');
model.sol('sol1').feature('s1').feature('se1').feature.create('ss2', 'SegregatedStep');
model.sol('sol1').feature('s1').feature('se1').feature.remove('ssDef');
model.sol('sol1').feature('s1').feature.remove('fcDef');

model.study('std1').feature('stat').set('initstudyhide', 'on');
model.study('std1').feature('stat').set('initsolhide', 'on');
model.study('std1').feature('stat').set('notstudyhide', 'on');
model.study('std1').feature('stat').set('notsolhide', 'on');

model.result.numerical.create('int1', 'IntVolume');
model.result.numerical.create('int2', 'IntVolume');
model.result.numerical('int1').selection.named('geom1_csel2_dom');
model.result.numerical('int1').set('probetag', 'none');
model.result.numerical('int2').selection.named('geom1_csel3_dom');
model.result.numerical('int2').set('probetag', 'none');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').name('Compile Equations: Stationary');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').feature('s1').set('control', 'stat');
model.sol('sol1').feature('s1').feature('se1').feature('ss1').set('segvar', {'comp1_s'});
model.sol('sol1').feature('s1').feature('se1').feature('ss2').set('segvar', {'comp1_i'});
model.sol('sol1').runAll;

model.result.numerical('int1').set('descr', 's*alpha1');
model.result.numerical('int1').set('expr', 's*alpha1');
model.result.numerical('int2').set('descr', 'i*alpha2');
model.result.numerical('int2').set('expr', 'i*alpha2');
model.result.numerical('int1').setResult;
model.result.numerical('int2').setResult;
J(1)=model.result.numerical('int1').getReal;
J(2)=model.result.numerical('int2').getReal;
end


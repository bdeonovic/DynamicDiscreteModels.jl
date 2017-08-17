
# m[x,y,x',y']=a[x,x']* b[x',y']
function hmm2ddm!(model::DynamicDiscreteModel, a, b)
	dx, dy = size(b)
	for jy in 1:dy
		for jx in 1:dx
			for iy in 1:dy
				for ix in 1:dx
					model.m[ix,iy,jx,jy] = a[ix,jx] * b[jx,jy]
				end
			end
		end
	end
end

#multiple dispatch will detech jacobian or no jacobian
function hmm2ddm!(model::DynamicDiscreteModel, a, b, ajac, bjac)
	dx,dy,dtheta=size(bjac)
	for jy in 1:dy
		for jx in 1:dx
			for iy in 1:dy
				for ix in 1:dx
					model.m[ix,iy,jx,jy] = a[ix,jx] * b[jx,jy]
					for itheta in 1:dtheta
						model.mjac[ix,iy,jx,jy,itheta] = ajac[ix,jx,itheta] * b[jx,jy] + 
                                             a[ix,jx] * bjac[jx,jy,itheta]
					end
				end
			end
		end
	end
end

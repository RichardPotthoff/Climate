import marimo

__generated_with = "0.15.0"
app = marimo.App()


@app.cell
def _():
    #load missing required files from GitHub, so the Notebook can be run in Google colab
    def wget(url,local_path='./'):
      filename=url.rsplit('/',1)[-1]
      import os
      if not os.path.exists(local_path+filename):
        import urllib3
        import certifi
        http=urllib3.PoolManager( cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
        with http.request('GET',url,preload_content=False) as r:
          with open(local_path+filename,'wb') as f:
            print(f'Downloading file: {url}')
            print(f'Saving file: {local_path}{filename}')
            f.write(r.read())
      else:
        print(f'"{local_path}{filename}" already exists. No need to download.')
    for filename in ('mollweide.py','fqs.py'):
      wget('https://raw.githubusercontent.com/RichardPotthoff/Climate/main/'+filename)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Constants""")
    return


@app.cell
def _():
    import numpy as np
    from numpy import pi,sin,cos, arccos as acos, arcsin as asin
    from matplotlib import pyplot as plt , cm as cmap 
    au=149_597_870_700 # [m] (by definition) astronomical unit: distance sun<->earth
    r_sun=695_700_000 # [m] radius of the sun
    T_sun=5_772.0 # [K] surface temperature of the sun
    r_earth=6_378_100 # [m] radius of the earth at equator
    g_earth=9.80665 # [m/s^2] acceleration due to gravity on earth
    cp_air=1005 # [J/(kg K)] specific heat capacity at constant pressure for dry air @15°C
    deg=pi/180 # 1 angle degree in radians (conversion factor)
    atm=101325 # [Pa] atmospheric pressure
    day=24*3600 # [s] length if day in seconds
    σ_SB=5.67037442e-8 # [J/(m^2 K^4)] Stefan-Boltzmann constant
    print(f'{r_sun=:7.3e}m, {T_sun=}K, {au=:7.3e}m')
    return (
        T_sun,
        acos,
        asin,
        atm,
        au,
        cmap,
        cos,
        cp_air,
        day,
        deg,
        g_earth,
        np,
        pi,
        plt,
        r_sun,
        sin,
        σ_SB,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # Black Body
    ## Solar irradiance
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""$$\dot q_{\rm au}=\sigma_\text{S-B}\cdot T_{\rm sun}^4\cdot\left(\frac{r_{\rm sun}}{1\cdot{\rm au}}\right)^2$$""")
    return


@app.cell
def _(T_sun, au, r_sun, σ_SB):
    q_au=σ_SB*T_sun**4*(r_sun/(1*au))**2
    print(f'{q_au=:0.0f}W/m^2')
    return (q_au,)


@app.cell
def _(atm, cp_air, g_earth):
    cp_atmosphere=cp_air*1*atm/g_earth
    print(f'{cp_atmosphere=:0.3e}J/(m^2 K)')
    return (cp_atmosphere,)


@app.cell
def _(cos, np, pi, plt):
    def Qsun(t):
      i,r=np.divmod(t+6,24)
      return i + (r>12)* 0.5*(1+cos(r*(pi/12)))
    tod=np.linspace(0,48,201)
    solrad=np.maximum(cos((tod-12)/24*2*pi),0)
    plt.plot(tod,solrad-1/pi,label="net heat flow: sun - IR")
    plt.plot(tod,Qsun(tod),zorder=20,label="cumulative solar heating")
    plt.plot((0,48),(0,-2),label="cumulative IR heat loss")
    plt.plot(tod,Qsun(tod)-tod/24,label="cumulative heat balance: solar in - IR out")
    plt.plot((0,48),(0,0),'k',zorder=10)
    plt.grid()
    plt.xlim((0,48))
    plt.xticks(range(0,49,3))
    plt.xlabel("time of day [h]")
    plt.ylabel("heat (total per day = 1)")
    plt.legend()
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Estimated daily temperature variation based on solar irradiance and heat capacity of atmosphere""")
    return


@app.cell
def _(acos, cos, cp_atmosphere, day, deg, pi, q_au, sin):
    def dT_dt(phi):
        return q_au * cos(phi) / (cp_atmosphere * pi)
    af = sin(acos(1 / pi)) - acos(1 / pi) * 1 / pi
    print(f'af={af!r}')
    for lat in (0, 30, 40, 50, 60, 70):
        print(f'lat={lat:2.0f}: dT_dt(lat*deg)*day = {dT_dt(lat * deg) * day:0.2f}K/day (ΔT={dT_dt(lat * deg) * day * af:0.2f}K with simultaneous cooling)')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Temperature of black-body and grey-body sphere at 1 au distance from sun""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""$$T_{\rm BB} = T_{\rm sun}\cdot\sqrt[\raisebox{-1pt}{$^4$}]{\frac{r_{\rm sun}^2}{{(1 \cdot \rm{au})}^2}\cdot\frac{\pi\cdot r_{\rm earth}^2}{4\cdot\pi\cdot r_{\rm earth}^2}}=T_{\rm sun}\cdot\sqrt[\raisebox{-1pt}{$^4$}]{\frac{r_{\rm sun}^2}{{(1 \cdot \rm{au})}^2}\cdot\frac{1}{4}}$$""")
    return


@app.cell
def _():
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""$$T_{\rm GB} = T_{\rm sun}\cdot\sqrt[\raisebox{-1pt}{$^4$}]{\frac{r_{\rm sun}^2}{{(1 \cdot \rm{au})}^2}\cdot\frac{(1-\alpha_{\rm cloud})\cdot\pi\cdot r_{\rm earth}^2}{(1-\alpha_{\rm cloud})\cdot4\cdot\pi\cdot r_{\rm earth}^2}}=T_{\rm sun}\cdot\sqrt[\raisebox{-1pt}{$^4$}]{\frac{r_{\rm sun}^2}{{(1 \cdot \rm{au})}^2}\cdot\frac{1}{4}}$$""")
    return


@app.cell
def _(T_sun, au, r_sun):
    T_BB = T_GB =T_sun * (r_sun**2/(1*au)**2 * (1/4))**(1/4)
    print(f'T_BB = T_GB = {T_BB:0.2f}K = {T_BB-273.15:0.2f}°C')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # Rotating planet (sun over equator)
    ## Temperature as function of latitude
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""$$T_{\mathrm rot}(\phi) = T_{\rm sun}\cdot \hspace{3pt} ^4 \hspace{-4pt} \sqrt{\frac{r_\mathrm{sun}^2}{{(1 \cdot \mathrm{au})}^2}\cdot\frac{{\cos}(\phi)}{\pi}}$$""")
    return


@app.cell
def _(mo):
    mo.md(r"""$$T_{\mathrm{rot}}(\phi) = T_{\mathrm{sun}} \cdot \raisebox{2pt}{$\hspace{3pt}^{4}$} \hskip -4pt \sqrt{\frac{r_{\mathrm{sun}}^2}{{(1 \cdot \mathrm{au})}^2} \cdot \frac{\cos(\phi)}{\pi}}$$""")
    return


@app.cell
def _(T_sun, au, cos, pi, r_sun):
    def T_rot(phi):
        return T_sun * (r_sun ** 2 / (1 * au) ** 2 * (cos(phi) / pi)) ** (1 / 4)
    return (T_rot,)


@app.cell
def _(cmap, cos, np, pi, sin):
    from mollweide import x_y2lat_lon, lat_lon2x_y

    def color_bar(ax):
        ax.set_xticks([])
        ax.yaxis.tick_right()
        ax.set_yticks([0, 1 / 3, 2 / 3, 1], labels=['-50°C', '  0°C', ' 50°C', '100°C'])
        ax.contourf(np.ones((1, 2)) * [[-50 + 273.15], [100 + 273.15]], cmap=cmap.jet, levels=100, vmin=220, vmax=220 + 150)

    def circle(r=1, Rx=None, Ry=None, tol=0.2):
        if Rx == None:
            Rx = r
        if Ry == None:
            Ry = r
        n = max(5, int(pi / (2 * tol / max(Rx, Ry)) ** 0.5) + 1)
        _phi = np.linspace(-pi, pi, n + 1)
        return (Rx * cos(_phi), Ry * sin(_phi))
    _R = 140

    def moll(ax, T_phi):
        ax.set_xticks([2 * _R * lon / 180 for lon in range(-180, 181, 30)], labels=[f'{lon:0.0f}' for lon in range(-180, 181, 30)])
        ax.set_yticks([_R * lat_lon2x_y(lat / 180 * pi, 0)[1] for lat in range(-90, 91, 30)], labels=[f'{lat:0.0f}' for lat in range(-90, 91, 30)])
        for lon in range(-180, 1, 30):
            ax.plot(*circle(_R, Rx=2 * _R * lon / 180, tol=0.2), 'k', zorder=20)
        ax.plot(*circle(1.01 * _R, Rx=1.01 * 2 * _R, tol=0.2), 'white', lw=4, zorder=15)
        for lat in range(0, 90, 30):
            x, y = lat_lon2x_y(lat / 180 * pi, pi, R=_R)
            ax.plot((-x, x), (y, y), 'k', zorder=20)
            ax.plot((-x, x), (-y, -y), 'k', zorder=20)
        for y in range(-_R, _R, 2):
            lat, lon = x_y2lat_lon(0, y, R=_R)
            x, _ = lat_lon2x_y(lat, -pi, R=_R)
            ax.plot((-x, x), (y, y), color=cmap.jet((T_phi(lat) - 220) / 150), lw=2.6, zorder=10)
        ax.set_title(f'Mollweide projection for "{T_phi.__name__}"')
    return color_bar, lat_lon2x_y, moll, x_y2lat_lon


@app.cell
def _(T_rot, deg):
    for _phi in reversed([0, 10, 20, 30, 40, 50, 60]):
        print(f'T_rot({_phi:2.0f}°Lat) = {T_rot(_phi * deg):0.2f}K, ={T_rot(_phi * deg) - 273.15:7.2f}°C')
    return


@app.cell
def _(T_rot, color_bar, moll, plt):
    _fig, (_ax, _cb) = plt.subplots(1, 2, figsize=(10, 17.5), gridspec_kw=dict(wspace=0.05, left=0, right=1, top=0.7, width_ratios=[2, 0.15]))
    _ax.set_aspect('equal')
    moll(_ax, T_rot)
    _cb.set_aspect(5.0)
    color_bar(_cb)
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Average temperature for a rotating planet""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""$$\overline{T_{\mathrm rot}}=T_{\mathrm sun}\cdot\sqrt[\raisebox{-1pt}{$^4$}]{\frac{r_\mathrm{sun}^2}{{(1 \cdot \mathrm{au})}^2}}\cdot\frac{\int_{-\frac \pi 2}^\frac \pi 2 \frac{\cos(\phi)^\frac 1 4}{\pi} \cdot \cos(\phi) \enspace {\mathrm d} \phi}{\int_{-\frac\pi 2}^\frac{\pi}{2} \cos(\phi)\enspace{\mathrm d}\phi}=T_{\mathrm sun}\cdot\sqrt[\raisebox{-1pt}{$^4$}]{\frac{r_\mathrm{sun}^2}{{(1 \cdot \mathrm{au})}^2}}\cdot{-\frac{4\cdot\sqrt[\raisebox{-1pt}{$^4$}]\pi\cdot\Gamma(\frac 1 8)}{15\cdot\Gamma({-\frac 3 8})}}$$""")
    return


@app.cell
def _(T_sun, au, pi, r_sun):
    from math import gamma
    T_rot_avg=T_sun * (r_sun**2/(1*au)**2)**(1/4)*-4/15*pi**(1/4)*gamma(1/8)/gamma(-3/8)
    print(f'T_rot_avg = {T_rot_avg:0.2f}K = {T_rot_avg-273.15:0.2f}°C')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # Tidally locked planet (or sun over pole)
    ## Temperature as function of latitude for a tidally locked planet
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""$$T_{\mathrm pol}(\phi) = T_{\mathrm sun}\cdot\sqrt[\raisebox{-1pt}{$^4$}]{\frac{r_\mathrm{sun}^2}{{(1 \cdot \rm{au})}^2}\cdot\frac{\sin(\phi)}{1}}$$""")
    return


@app.cell
def _(T_sun, au, deg, r_sun, sin):
    def T_pol(phi):
        return 0 if phi < 0 else T_sun * (r_sun ** 2 / (1 * au) ** 2 * (sin(phi) / 1)) ** (1 / 4)
    for _phi in [90, 80, 70, 60, 50, 40, 30, 20, 10, 0, -10]:
        print(f'T_rot({_phi:3.0f}°Lat) ={T_pol(_phi * deg):7.2f}K ={T_pol(_phi * deg) - 273.15:7.2f}°C')
    return (T_pol,)


@app.cell
def _(T_pol, color_bar, moll, plt):
    _fig, (_ax, _cb) = plt.subplots(1, 2, figsize=(10, 17.5), gridspec_kw=dict(wspace=0.05, left=0, right=1, top=0.7, width_ratios=[2, 0.15]))
    _ax.set_aspect('equal')
    moll(_ax, T_pol)
    _cb.set_aspect(5.0)
    color_bar(_cb)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Average temperature for a tidally locked planet""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""$$\overline{T_{\rm pol}}=T_{\mathrm sun}\cdot\sqrt[\raisebox{-1pt}{$^4$}]{\frac{r_\mathrm{sun}^2}{{(1 \cdot \rm{au})}^2}}\cdot\frac{\int_0^\frac \pi 2 \frac{\sin(\phi)^\frac 1 4}{1} \cdot \cos(\phi) \enspace {\rm d} \phi}{\int_{-\frac\pi 2}^\frac{\pi}{2} \cos(\phi)\enspace{\mathrm d}\phi}=T_{\mathrm sun}\cdot\sqrt[\raisebox{-1pt}{$^4$}]{\frac{r_\mathrm{sun}^2}{{(1 \cdot \mathrm{au})}^2}}\cdot\frac 2 5$$""")
    return


@app.cell
def _(T_sun, au, r_sun):
    T_pol_avg=T_sun * (r_sun**2/(1*au)**2)**(1/4)*2/5
    print(f'T_pol_avg = {T_pol_avg:0.2f}K = {T_pol_avg-273.15:0.2f}°C')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # Appendix
    ## Side-by-side comparison of temperature maps for a tidally locked and a rotating planet
    """
    )
    return


@app.cell
def _(T_pol, T_rot, asin, cmap, color_bar, cos, deg, np, pi, plt, sin):
    def circle_1(r=1, Rx=None, Ry=None, tol=0.2):
        if Rx == None:
            Rx = r
        if Ry == None:
            Ry = r
        n = max(5, int(pi / (2 * tol / max(Rx, Ry)) ** 0.5) + 1)
        _phi = np.linspace(-pi, pi, n + 1)
        return (Rx * cos(_phi), Ry * sin(_phi))

    def flat_earth(ax, T_phi):
        ax.set_aspect('equal')
        ax.set_axis_off()
        R = 140
        ax.set_xlim(-R * 1.01, R * 1.01)
        ax.set_ylim(-R * 1.01, R * 1.01)
        for lon in range(0, 180, 30):
            ax.plot((-R * cos(lon * deg), R * cos(lon * deg)), (-R * sin(lon * deg), R * sin(lon * deg)), 'k', zorder=10)
        for lat30 in range(-90, 89, 30):
            ax.plot(*circle_1(R / 2 * (2 - 2 * sin(lat30 * deg)) ** 0.5, tol=0.2), 'r' if lat30 == 0 else 'k' if T_phi(lat30 * deg) > 100 or lat30 == -90 else 'gray', zorder=10)
        for r in range(2, R, 2):
            _phi = asin(1 - (r / (R / 2 ** 0.5)) ** 2)
            ax.plot(*circle_1(r, tol=0.2), color=cmap.jet((T_phi(_phi) - 220) / 150), lw=2.6, zorder=5)
        ax.set_title(f'Equal-Area Temperature Map for "{T_phi.__name__}"')
    plt.close()
    _fig, _axes = plt.subplots(1, 3, figsize=(10.0, 9.0), gridspec_kw=dict(wspace=0.05, left=0, right=1, top=0.7, width_ratios=[1, 1, 0.15]))
    flat_earth(_axes[0], T_pol)
    flat_earth(_axes[1], T_rot)
    _axes[2].set_aspect(5.0)
    color_bar(_axes[2])
    _fig
    return (circle_1,)


@app.cell
def _(
    T_pol,
    T_rot,
    circle_1,
    cmap,
    color_bar,
    lat_lon2x_y,
    pi,
    plt,
    x_y2lat_lon,
):
    _R = 140

    def moll_1(ax, T_phi):
        ax.set_xticks([2 * _R * lon / 180 for lon in range(-180, 181, 30)], labels=[f'{lon:0.0f}' for lon in range(-180, 181, 30)])
        ax.set_yticks([_R * lat_lon2x_y(lat / 180 * pi, 0)[1] for lat in range(-90, 91, 30)], labels=[f'{lat:0.0f}' for lat in range(-90, 91, 30)])
        for lon in range(-180, 1, 30):
            ax.plot(*circle_1(_R, Rx=2 * _R * lon / 180, tol=0.2), 'k', zorder=20)
        ax.plot(*circle_1(1.01 * _R, Rx=1.01 * 2 * _R, tol=0.2), 'white', lw=4, zorder=15)
        for lat in range(0, 90, 30):
            x, y = lat_lon2x_y(lat / 180 * pi, pi, R=_R)
            ax.plot((-x, x), (y, y), 'k', zorder=20)
            ax.plot((-x, x), (-y, -y), 'k', zorder=20)
        for y in range(-_R, _R, 2):
            lat, lon = x_y2lat_lon(0, y, R=_R)
            x, _ = lat_lon2x_y(lat, -0.99 * pi, R=_R)
            ax.plot((-x, x), (y, y), color=cmap.jet((T_phi(lat) - 220) / 150), lw=2.6, zorder=10)
        ax.set_title(f'Mollweide projection for "{T_phi.__name__}"')
    _fig, ((_ax0, _cb0), (_ax1, _cb1)) = plt.subplots(2, 2, figsize=(10, 17.5), gridspec_kw=dict(wspace=0.05, left=0, right=1, top=0.7, width_ratios=[2, 0.15]))
    _ax0.set_aspect('equal')
    moll_1(_ax0, T_rot)
    _ax1.set_aspect('equal')
    moll_1(_ax1, T_pol)
    _cb0.set_aspect(5.0)
    color_bar(_cb0)
    _cb1.set_aspect(5.0)
    color_bar(_cb1)
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## NASA image of earth""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""![IMG_0097.AVIF](https://raw.githubusercontent.com/RichardPotthoff/Climate/main//IMG_0097.AVIF)""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Image of the Chimborazo""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""![Chimborazo](https://raw.githubusercontent.com/RichardPotthoff/Climate/main//IMG_0099.WEBP)""")
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()

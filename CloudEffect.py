# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "marimo",
#     "matplotlib==3.10.6",
#     "numpy==2.2.6",
# ]
# ///

import marimo

__generated_with = "0.16.5"
app = marimo.App(width="medium")

with app.setup(hide_code=True):
    # Initialization code that runs before all other cells
    import marimo as mo
    from matplotlib import pyplot as plt
    import numpy as np
    from collections import namedtuple
    from numpy import sin,cos,pi
    import re
    from typing import Self
    import sys
    au=149_597_870_700 # [m] (by definition) astronomical unit: distance sun<->earth
    r_sun=695_700_000 # [m] radius of the sun
    T_sun=5_772.0 # [K] surface temperature of the sun
    r_earth=6_378_100 # [m] radius of the earth at equator
    g_earth=9.80665 # [m/s^2] acceleration due to gravity on earth
    cp_air=1005 # [J/(kg K)] specific heat capacity at constant pressure for dry air @15°C
    deg=pi/180 # 1 angle degree in radians (conversion factor)
    atm=101325 # [Pa] atmospheric pressure
    day=24*3600 # [s] length if day in seconds
    σ_SB=5.67037442e-8 # W/(m**2 K)

    Planet_csv="""
    name,R_orbit,r_per,r_ap
    str,float,float,float
    BB,0,0,0
    Mercury, 57909227000.0, 46001009000.0, 69817445000.0
    Venus, 108209475000.0, 107476170000.0, 108942780000.0
    Earth, 149598262000.0, 147098291000.0, 152098233000.0
    Mars, 227943823500.0, 206655215000.0, 249232432000.0
    Jupiter, 778340821000.0, 740679835000.0, 816001807000.0
    Saturn, 1426666422000.0, 1349823615000.0, 1503509229000.0
    Uranus, 2870658186000.0, 2734998229000.0, 3006318143000.0
    Neptune, 4498396441000.0, 4459753056000.0, 4537039826000.0
    """ 

    def csv2tuples(csv,tuplename="NT"):
        tm={"str":str ,"float":float}
        records=([field.strip() for field in line.split(',')] for line in csv.strip().splitlines())
        names=next(records)
        types=[tm[typestr] for typestr in next(records)]
        NT=namedtuple(tuplename,names)
        return {types[0](record[0]):NT(*(T(value) for T,value in zip(types,record))) for record in records}

    Planets=csv2tuples(Planet_csv)

    Length_units={
        "m":1.0, # m
       "km":1000.0, # m
       "au":au,
    }
    units={}
    units|=Length_units

    def sci_tex(value: float, fmt: str = "0.3f", unit: str = None, style: str = None, ) -> str:

        from math import floor,log10
        if style==None:
            style="s"#default if scientific E format  
        style=style.lower()
        if style in "s31" and ("e" in fmt or "E" in fmt):
            fmt=fmt.replace("e","f").replace("E","F")
            if value == 0:
                tex=rf"{value:{fmt}}"
            else:
                abs_value = abs(value)
                logval=log10(abs_value)
                exp=floor(logval)
                mantlog=logval-exp
                expstep=3 if style =="3" else 1
                expsteps,mantcorr=divmod(exp,expstep)
                exp=expsteps*expstep
                mantlog=mantlog+mantcorr
                mantissa=10**mantlog
                if value<0:
                    mantissa*=-1
                tex=rf"{{{mantissa:{fmt}} \times 10^\mathrm{{{exp}}}}}"
        else:
            tex=rf"\text{{{value:{fmt}}}}"

        if unit:
            processed_unit = (unit
                .replace("**", "^")
                .replace("*", r"\; ")
                .replace("°C",r"{^{\boldsymbol{\circ}}}\mathrm{C}")
                )
            enu_deno = processed_unit.split(':', 1)
            if len(enu_deno)==1:
                tex += rf" \,\mathrm{{{enu_deno[0]}}}"
            else:
                tex += rf" \,\mathrm{{\frac{{{enu_deno[0]}}}{{{enu_deno[1]}}}}}"

        return tex

    class TeXfloat(float):
        """
        Subclass of float with custom 'T' format specifier for LaTeX scientific/engineering notation.

        Usage:
            >>> tf = TeXFloat(1234.567)
            >>> f"{tf:.3engT[m]}"  # '{1.235 \\times 10^{3}} \\,\\mathrm{m}'
            >>> f"{tf:.2eT[kg*m/s**2]}"  # '{1.23 \\times 10^{3}} \\,\\mathrm{\\frac{kg \\cdot m}{s^{2}}}'
            >>> f"{tf:.3T[]}"  # '{1.235 \\times 10^{3}}' (no unit)
            >>> tf * 2  # 2469.134 (plain float ops)
            >>> f"{tf + 100:.0T[m/s]}"  # '{1 \\times 10^{3}} \\,\\mathrm{\\frac{m}{s}}' (prec=0, default eng)

        Spec format: :.prec[eng|e]T[unit]
        - prec: Decimal places (default 3).
        - eng/e: Mode before T (default eng).
        - unit: Bracketed string (default None); supports /, **, *.
        """
        def __radd__(self, other: float) -> Self:   return(self.__class__(super().__radd__(other)))
        def __rsub__(self, other: float) -> Self:   return(self.__class__(super().__rsub__(other)))
        def __rmul__(self, other: float) -> Self:   return(self.__class__(super().__rmul__(other)))
        def __rtruediv__(self, other: float)->Self: return(self.__class__(super().__rtruediv__(other))) 
        def __format__(self, spec: str) -> str:
            pattern =  r'^(?P<format>[^Tt]*)(?:(?P<sep>[tT])(?P<style>[^[])?(?:\[(?P<unit>.*))\])?$'
            fmt, unit, style=spec, None, None
            match=re.match(pattern,spec)
            if match:      
                gd =match.groupdict()
                if gd["sep"]:
                    fmt=gd["format"]
                    style=gd["style"]
                    unit=gd["unit"]
            # Fallback to standard float format
            return sci_tex(float(self), fmt=fmt, style=style, unit=unit)

    #constants that can be used for converting float to TeXfloat
    TeX_0=TeXfloat(0.0)
    TeX_1=TeXfloat(1.0)
    TeX_C=TeXfloat(273.15)
    TeX_p=TeXfloat(1e-12) 
    TeX_n=TeXfloat(1e-9)
    TeX_μ=TeXfloat(1e-6)
    TeX_m=TeXfloat(1e-3)
    TeX_c=TeXfloat(1e-2)
    TeX_h=TeXfloat(1e2)
    TeX_k=TeXfloat(1e3)
    TeX_M=TeXfloat(1e6)
    TeX_G=TeXfloat(1e9)
    TeX_T=TeXfloat(1e12)


@app.cell
def _():
    mo.md(
        f"""
    # Table of Contents
    {mo.outline()}
    """
    )
    return


@app.cell(hide_code=True)
def _(TeX_Macro):
    mo.md(f"{TeX_Macro}"r"""
    # Symbols and Constants
    | Symbol | Definition |
    | ---:|:--- |
    | $\rmsub{T}{BB}$ | Surface Temperature of a Black Body ${}$ |
    | $\rmsub{T}{sun}$ | Surface temperature of the sun ( """f"{ T_sun:0,.0f}"r""" K ) ${}$ |
    | $\rmsub{r}{sun}$ | Radius of the sun ( """f"{r_sun/1000:0,.0f} km )"r""" ${}$ |
    | $\rmsub{R}{orbit}$ | Distance from the sun ${}$ |
    | $\rmsub{g}{\,Eearth}$ | Acceleration due to gravity on earth ( """f"{g_earth:0.2f}"r"""$\,\mathrm{m} / \mathrm{s}^2$ )|
    | $\phi$ | Latitude ${}$ |
    | $\mathrm{au}$ | Astronomical unit ( distance sun$\,\leftrightarrow\,$earth = """f"{au/1000:0,.0f} km )"r"""
    | $\omega$|cloud cover (fraction) ${}$ |
    | $\Gamma$|Temperature lapse rate ${}$ |
    | $\rmsub{\sigma}{SB}$| The Stefan-Boltzmann constant ( """f"${σ_SB+TeX_0:0.3eT3[W/(m**2*K)]}$"r""" ) |
    """ )
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""# Earth as a Black Body""")
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    ## Uniform Temperature

    The temperature of a black body in orbit around the sun depends on the temperature of the sun, the size of the sun, and distance of the planet from the sun: 

    $$T_\mathrm{BB}=T_\mathrm{sun}\cdot\sqrt[\raisebox{-1pt}{$^4$}]{\frac{r_\mathrm{sun}^2}  {R_\mathrm{orbit}^{\,2}}\cdot\frac{1}{4}}$$

    The factor $1/4$ is the ratio of the projected area of earth (its shaddow area), divided by the surface area (which radiates infrared radiation into space).
    """
    )
    return


@app.cell(hide_code=True)
def _():
    #States
    get_R_orbit,set_R_orbit=mo.state(au)
    return get_R_orbit, set_R_orbit


@app.cell(hide_code=True)
def _(set_R_orbit):
    #Selectors
    def update_R_orbit(*args,**kwargs):
        new_R_orbit=Planet_selector.value[Orbit_selector.value]
        if new_R_orbit != 0.0:
            set_R_orbit(new_R_orbit)

    Planet_selector=mo.ui.dropdown(options=Planets, value="Earth", on_change=update_R_orbit) 
    Orbit_options={key:Planets["Earth"]._fields.index(value) for key,value in {"orbit avg.":"R_orbit","perihelion":"r_per","aphelion":"r_ap"}.items()}
    Orbit_selector=mo.ui.dropdown(options=Orbit_options,value="orbit avg.", on_change=update_R_orbit)
    R_orbit_unit=mo.ui.dropdown(options=Length_units,value="au")
    return Orbit_selector, Planet_selector, R_orbit_unit


@app.cell(hide_code=True)
def _(R_orbit_unit, get_R_orbit, set_R_orbit):
    #Edit number
    R_orbit_number=mo.ui.number(value=get_R_orbit()/R_orbit_unit.value, on_change=lambda R_au:set_R_orbit(R_au*R_orbit_unit.value))
    return (R_orbit_number,)


@app.cell(hide_code=True)
def _(get_R_orbit):
    #Calculations 
    T_BB=T_sun*(r_sun**2/get_R_orbit()**2*1/4)**(1/4)
    return (T_BB,)


@app.cell(hide_code=True)
def _(
    Orbit_selector,
    Planet_selector,
    R_orbit_number,
    R_orbit_unit,
    T_BB,
    get_R_orbit,
):
    mo.md(
        rf"""
    {Planet_selector}   $\;R\;${Orbit_selector}$\;=\;${R_orbit_number if not Planet_selector.value[Orbit_selector.value] else get_R_orbit()/R_orbit_unit.value:0.6} {R_orbit_unit}$\,=\,{get_R_orbit()*TeX_1:0.3eT3[m]}$

    $$T_\mathrm{{ {Planet_selector.selected_key} }}
    ={T_sun+TeX_0:.3eT3[K]} \cdot \sqrt[\raisebox{{-1pt}}{{$^4$}}]{{\frac{{({r_sun+TeX_0:,.3eT3[m]})^{{\,2}}}} {{({get_R_orbit()+TeX_0:0,.3eT3[m]})^{{\,2}}}}\cdot\frac{{1}}{{4}}}}
    ={T_BB+TeX_0:0.1fT[K]}
    ={T_BB-TeX_C:0.1fT[°C]}$$
    """
    )
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""$$q_\mathrm{{orbit}}=\sigma_\mathrm{{SB}} \cdot T_{{sun}}^4\cdot \frac{{r_\mathrm{{sun}}^{{\,2}}}} {{R_\mathrm{{orbit}}^{{\,2}}}}$$""")
    return


@app.cell(hide_code=True)
def _(
    Orbit_selector,
    Planet_selector,
    R_orbit_number,
    R_orbit_unit,
    get_R_orbit,
):
    mo.md(
        rf"""
    {Planet_selector}   $\;R\;${Orbit_selector}$\;=\;${R_orbit_number if not Planet_selector.value[Orbit_selector.value] else get_R_orbit()/R_orbit_unit.value:0.6} {R_orbit_unit}$\;=\,{sci_tex(get_R_orbit(),"0.3E","m","3")}$

    $$
    q_\mathrm{{\,{Planet_selector.selected_key}}} = {sci_tex(σ_SB,"0.3e","W:m^2*K**4","3")} \cdot ({sci_tex(T_sun,"0.3e","K")})^4
    \cdot \frac{{({sci_tex(r_sun,"0.3e","m","3")})^{{\,2}}}} {{({sci_tex(get_R_orbit(),"0.3e","m","3")})^{{\,2}}}}
    = {sci_tex(

    σ_SB * T_sun**4 * r_sun**2 / get_R_orbit()**2 ,

    "0.3e","W:m^2","3")}
    $$
    """
    )
    return


@app.cell
def _():
    mo.md(
        r"""
    ## Effect of Clouds on the Temperature of the Planet
    ### Grey Body
    The eart is not a black body, so a correction has to be made for the light that is reflected and not absorbed. The albedo of earth is for the most part caused by clouds that have a high reflectivity compared to the oceans and the land. The ice at the poles has also a high reflectivity for sunlight, but there is not much light that could be reflected at the poles. 
    To keep the model simple, I assume (for now) that the albedo of the earth is caused by the clouds. 
    Another correction has to be made for the infrared radiation from the surface that is blocked by the same clouds. 
    For the case that the fraction of daylight relected by clouds, $\  \omega$, is the same as the fraction of the infrared light that is reflected back to earth we can  calculate the temperature of this "grey body": 

    $$T_\mathrm{GB}=T_\mathrm{sun}\cdot\sqrt[\raisebox{-1pt}{$^4$}]{
    \frac{1-\omega}  {1-\omega}
    \cdot\frac{r_\mathrm{sun}^2}  {R_\mathrm{orbit}^{\,2}}
    \cdot\frac{1}{4}}$$

    It turns out, that the temperature of a grey body is the same as the temperature of a grey body, because the $({1-\omega})$ terms in the enumerator and denominator cancel out.
    To make the model more realistic, we can add the infrared radiation emitted by the clouds back in, but the amount of radiation emitted by the clouds depends on the temperature, and the temperature is a function of the altitude: High clouds will be colder and emit less infraded radiation than low clouds fog near the ground.
    """
    )
    return


app._unparsable_cell(
    r"""
    mo.md(
        f\"\"\"
    {TeX_Macro}\" r\\"\"\"
    $$T_\mathrm{GB}=T_\mathrm{sun}\cdot\sqrt[\raisebox{-1pt}{$^4$}]{
    \frac{1-\omega}  {1-\omega+\omega \cdot\epsilon}
    \cdot\frac{r_\mathrm{sun}^2}  {R_\mathrm{orbit}^{\,2}}
    \cdot\frac{1}{4}}$$

    $$ \epsilon
    =\frac{\rmsub{q}{cloud}}{\rmsub{q}{surf}}
    =\left(\frac{\rmsub{T}{cloud}}{\rmsub{T}{surf}}\right)^4
    $$

    $$ \rmsub{T}{cloud}=\rmsub{T}{surf} + \Gamma \cdot \rmsub{h}{cloud}$$

    $$ \epsilon
    =\left(\frac{\rmsub{T}{surf} + \Gamma \cdot \rmsub{h}{cloud}}{\rmsub{T}{surf}}\right)^4
    =\left(1+\frac{\Gamma}{\rmsub{T}{surf}} \cdot \rmsub{h}{cloud}\right)^4
    $$ 

    Since $\Gamma$ is negative, we could get negative argument values for very high values of $\rmsub{h}{cloud}$. By limiting the argument to values $\gt$ 0, we can limit the $\epsilon$ to a range between $\epsilon=1$ (\ for $\rmsub{h}{cloud}=0$\ ) and $\epsilon=0$ (\ for $\rmsub{h}{cloud} \gt \rmsub{T}{surf}/\;{\Gamma}$\ ):

    $$ \epsilon (\rmsub{h}{cloud})
    =\left[\mathrm{max}\left(0,\, 1+\frac{\Gamma}{\rmsub{T}{surf}} \cdot \rmsub{h}{cloud}\right)\right]^4
    $$ 

    For Earth the value of $\rmsub{T}{surf}/\;{\Gamma}$ is approximately $288.15\;\mathrm{K}\; /\; 6.5\mathrm{\frac{K}{km}}=44.33\;\mathrm{km}$. This is much higher than the highest tropospheric clouds. 
    \"\"
    \"\"\"
    )
    """,
    name="_"
)


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    ## Temperature as a Function of the Latitude

    Assuming that the Temperature is the same across the whole surface of the black body is not very realistic for large bodies like the earth: The poles are much colder than the equator because the rays from the sun hit the surface at a shallower angle.

    A more realistic assumption would be that the temperature is the same only for a given latitude, since the rotation of the planet helps to equalize the temperature if the planet rotates fast enough. 

    The factor $1/4$ can then be replaced by the ratio $\cos(\phi)/\pi$. This the projected area of a conical frustum divided by its mantle surface area:

    $$T_\mathrm{BB}(\phi) = T_\mathrm{sun}\cdot\sqrt[\raisebox{-1pt}{$^4$}]{\frac{r_\mathrm{sun}^2}  {R_\mathrm{orbit}^{\,2}}\cdot\frac{\cos(\phi)}{\pi}}$$

    To keep the equation simple, it was assumed that the axis of rotation is perpendicular to the direction of the sun. This is strictly true only at the spring and fall equinoxes.
    """
    )
    return


@app.cell
def _():
    mo.md(
        r"""
    ### Cloud Temperature as Function of Altitude
    The main reason for the vertical temperature gradient in the atmosphere is gravity, as can be deduced trom the equation for the temperature laps rate $\Gamma$ :


    $$ \Gamma = \frac{\mathrm dT}{\mathrm dh} = -\frac{g}{c_{p}}$$

    In this equation, $\ g$ is the acceleration due to gravity, and $c_p$ the specific heat capacity of the air at constant pressure. For dry air the temperature lapse rate can be calculated as follows:
    """
    )
    return


@app.cell
def _():
    mo.md(
        rf"""
    $$\Gamma_\mathrm{{dry\ air}}
    =-\frac{{{g_earth+TeX_0:.2fT[m/s^2]}}}{{{cp_air/TeX_k:0.3fT[kJ/kg*K]}}}
    =-{g_earth/cp_air*TeX_k:0.2fT[°C/km]}$$
    """
    )
    return


@app.cell
def _():
    mo.md(
        rf"""
    This is the calculated value for *dry* air. For humid air, the heat capacity of the air is larger, because of the latent heat of water vapor. This leads to a lower value for the temperature lapse rate. Thermal convection in the atmosphere also reduces the temperature difference between high altitudes and low altitudes, because air from warmer regions of the earth tend to layer on top of colder air from colder regions, thus further reducing the thermal gradient. A commonly used average value for the temperature lapse rate in earth's atmosphere is: 

    $$
    \Gamma 
    \approx -{TeXfloat(6.5):0.1fT[°C/km]}$$
    """
    )
    return


@app.cell(hide_code=True)
def _():
    h_Trop=12000 #m
    dT_dh=-6.5/1000 #K/m
    h=np.linspace(0,h_Trop,121)
    for T_ref in reversed((-20,0,20,40)):
       f_corr=(1+dT_dh/(T_ref+273.15)*h)**4
       plt.plot(h,f_corr,label=f"{T_ref=}°C")

    dT_dh/T_ref
    plt.legend()
    plt.gca()
    return


@app.cell
def _():
    return


@app.cell
def _():
    mo.md(r"""#Appendix""")
    return


@app.cell
def _():
    mo.md(r"""## TeX Macro Definitions""")
    return


app._unparsable_cell(
    r"""
    mo.md(
        f\"\"\"
    {TeX_Macro}\"
    r\\"\"\"
    $$\dd T z ,\ \Grad A ,\ \rot {\mathbf v},\ \rmsub{T}{cloud}$$
    \"\"
    \"\"\"
    )
    """,
    name="_"
)


@app.cell
def _():
    #This is a setup cell for global LaTeX definitions and css styles
    #This cell should run only once if it is not changed, because it does
    #not have any inputs (except mo.md).
    #Add new definitions to the top of this list as the re-definitions of
    #the alrealy existing commands will cause an error that prevents the
    #statements below this line from being executed 
    _TeX_Macro=r"""$$
    \newcommand{\rmsub}[2]{#1_\mathrm{#2}}
    \newcommand{\dd}[2]{\frac{\mathrm d\: #1}{\mathrm d\: #2}}
    \newcommand{\Grad}[1]{\mathrm{Grad} \large( #1 \large)}
    \newcommand{\rot}[1]{\mathrm{curl} \left( #1 \right)}
    $$
    """
    #changed markdown style
    _css="""
    <style>
    /* Minimal scientific table: horizontal lines only, no grids */
    .prose table {
      border-collapse: collapse;
      width: 100%;
      font-size: 14px;
      margin: 1em 0;
      border-top: 1px solid #ddd;
      border-bottom: 1px solid #ddd;
    }

    .prose table th,
    .prose table td {
      border: none; /* No verticals or cell borders */
      padding: 4px 6px;
      text-align: left;
      background-color: white;
      vertical-align: top;
    }

    .prose table th {
      border-bottom: 2px solid #333; /* Header underline */
      font-weight: bold;
    }

    /* Remove even-row striping */
    .prose table tbody tr:nth-child(even) {
      background-color: transparent;
    }

    /* Dark mode support */
    @media (prefers-color-scheme: dark) {
      .prose table {
        border-top-color: #4b5563;
        border-bottom-color: #4b5563;
      }
      .prose table th,
      .prose table td {
        background-color: #1f2937;
      }
      .prose table th {
        border-bottom-color: #d1d5db;
      }
    }
    </style>
    """
    #The global variable "TeX_Macro" is used to create dependencies between this cell
    #and the cells using the macros. This forces this cell to run before the cells that 
    #use the macros defined here.
    #The easiest way to creaate a dependency to this cell is to reference the "TeX_Maxro"
    #variable by using it in a markdown f string: mo.md(f"{TeX_Macro} ..."). The value of 
    #TeX_Macro is an empty string, so it will not show up in the rendered markdown text.
    TeX_Macro=""
    mo.md(_css+_TeX_Macro) #run the macro definitions by rendering the _TeX_Macro string.
    return (TeX_Macro,)


if __name__ == "__main__":
    app.run()

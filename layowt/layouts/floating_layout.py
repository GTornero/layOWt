"""
This module defines the FloatingLayout class.
This class can be used to design and represent floating turbines and their mooring lines.
"""

from typing import Callable

import fiona
import geopandas as gp
import numpy as np
import pyproj
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from shapely.geometry import MultiLineString, MultiPolygon, Polygon, mapping
from shapely.ops import unary_union

from layowt.grids import Grid
from layowt.layouts.layout import Layout, _is_epsg


class FloatingLayout(Layout):
    """
    This class is responsible for holding all of the data structures required \
    in order to define a floating project scenario and to generate the coordinate pairs \
    of the points that make up the grid. \
    
    This Layout subclass adds onto the Layout class by adding functions and methods to \
    generate valid mooring layout geometries. All attributes, methods, and properties of the Layout \
    class are available to the FloatingLayout class. \

    Multiple numerical, spatial and graphical methods are available to \
    manipulate the scenario o conform to constraint and set properties. \
    
    Parameters
    ----------
    grid : Grid | None, optional
        Grid object used to define the Layout, by default None.
    epsg : int | None, optional
        EPSG code for the layout coordinates, by default None.
    areas : list[Polygon  |  MultiPolygon] | None, optional
        List of Shapely Polygons or MultiPolygons used to define the lease
        areas, by default None.
    exclusions : list[Polygon  |  MultiPolygon] | None, optional
        List of Shapely Polygons or MultiPolygons used to define exclusion
        zones where turbines cannot be placed, by default None.
    
    See Also
    --------
    Polygon : Shapely Polygon object.
    MultiPolygon : Shapely MultiPolygon object.
    MultiLineString : Shapely MultiLineString object.
    
    Attributes
    ----------
    grid : Grid
        Attribute: Underlying grid to the layout. None for layouts not created from a grid.
    epsg : int
        Attribute: `EPSG <https://epsg.io/>`_ code of the layout coordinate system.
    area : list[Polygon | MultiPolygon]
        Attrbute: List of buildable area geometries of the layout.
    exclusion : list[Polygon | MultiPolygon]
        Attribute: List of exclusion area geometries of the layout.
    bathymetry_path : str
        Attribute: filepath of the layout's bathymetry dataset.
    bathymetry_limits: tuple(float, float)
        Attribute: Tuple of minimum and maximum allowable bathymetry values for valid layout positions.
    bathymetry_sign : str
        Attribute: Sign ("+" or "-") in which bathymetry values are represented.
    bathymetry_drop_na : bool
        Attributes: If layout positions at locations of missing bathymetry have been removed.
    bathymetry_epsg : int
        Attribute: EPSG code of the bathymetry raster file.
    bathymetry_band : int
        Attribute: Band number of the bathymetry raster file being used.
    samples : list[float]
        Attribute: List of sampled bathymetry values at the layout positions.
    geom : MultiPoint
        Attribute: MultiPoint geometry representation of the layout positions.
    mooring_geom : list[MultiLineString]
        Attribute: A list of MultiLineString geometry representations for the mooring lines at each
        turbine position.
    n_moorings : int
        Attribute: Number of mooring lines per turbine position.
    mooring_radius_func : float | Callable
        Attribute: Function used to calculate the mooring radius at each turbine position. Can be
        a constant float or a Callable that takes local turbine bathymetry as its sole argument.
    platform_bearing: float
        Attribute: Cardinal bearing of the floating platforms in degrees from True North.
    mooring_bearings: list[float]
        Attribute: List of cardinal bearings of the platform mooring lines in the global frame
        of reference.
    mooring_offsets: list[tuple(float, float)]
        Attribute: List of tuple coordinate pairs for the mooring fairlead locations with respect
        to the turbine location in the local platform reference frame.
    mooring_origins: list[list[tuple(float, float)]]
        Attribute: List of lists of tuple coordinate pairs. Each list has tuple coordinate pairs
        for the origins for the mooring fairleads for each turbine platform in the global frame of
        reference.
    mooring_anchors: list[list[tuple(float, float)]]
        Attribute: List of lists of tuple coordinate pairs. Each list has tuple coordinate pairs
        for the end points for the mooring anchors for each turbine platform in the global frame of
        reference.
    mooring_radii: np.array
        Attribute: Numpy array of the mooring radius for each turbine platform.
    mooring_area: Polygon | MultiPolygon
        Attribute: The polygon(s) geometry which must contain the mooring lines.
    mooring_exclusions: Polygon | MultiPolygon
        Attribute: The polygon(s) geometries which the mooring lines cannot fall into.
    mooring_buffers: MultiPolygon
        Attribute: The MultiPolygon geometries of the mooring buffers.
    """
    def __init__(self,
        grid: Grid | None = None,
        epsg: int | None = None,
        areas: list[Polygon | MultiPolygon] | None = None,
        exclusions: list[Polygon | MultiPolygon] | None = None
    ):
        super().__init__(grid=grid, epsg=epsg, areas=areas, exclusions=exclusions)
        self.mooring_geom = None
        self.n_moorings = None
        self.mooring_radius_func = None
        self.platform_bearing = None
        self.mooring_bearings = None
        self.mooring_offsets = None
        self.mooring_origins = None
        self.mooring_anchors = None
        self.mooring_radii = None
        self.mooring_area = None
        self.mooring_exclusions = None
        self.mooring_buffers = None
    
    def _set_platform_bearing(self, platform_bearing: float | None) -> None:
        """
        Private setter method that sets the platform_bearing instance attribute.

        Parameters
        ----------
        platform_bearing : float | None
            Cardinal bearing of the platform from True North in the global reference frame.
        """
        if platform_bearing is None:
            self.platform_bearing = 0
        else:
            self.platform_bearing = platform_bearing % 360
    
    def _set_mooring_offsets(self, mooring_offsets: list[tuple[float, float]] | None) -> None:
        """
        Private method to calculate the mooring offsets after the rotating the platform by the \
        platform_bearing.

        Parameters
        ----------
        mooring_offsets : list[tuple[float, float]] | None
            List of coordinate pairs of the moorings in the platform local reference frame before
            applying the platform_bearing rotation.

        Raises
        ------
        ValueError
            The number of moorings does not match the number of coordinate pairs passed.
        """
        if mooring_offsets is None:
            self.mooring_offsets = self._rotate_mooring_offests([(0, 0)] * self.n_moorings)
        elif len(mooring_offsets) == self.n_moorings:
            if self.platform_bearing == 0:
                self.mooring_offsets = mooring_offsets
            else: self.mooring_offsets = self._rotate_mooring_offests(mooring_offsets)
        else:
            msg = f"The length of mooring_offsets must be the same as n_moorings. Has {len(mooring_offsets)} instead of the expected {self.n_moorings}."
            raise ValueError(msg)
    
    def _rotate_mooring_offests(self, mooring_offsets: list[tuple[float, float]]) -> list[tuple[float, float]]:
        """
        Private method to rotate the mooring fairlead coordinates.

        Parameters
        ----------
        mooring_offsets : list[tuple[float, float]]
            List of coordinate pairs of the moorings in the platform local reference frame before
            applying the platform_bearing rotation.

        Returns
        -------
        list[tuple[float, float]]
             List of coordinate pairs of the moorings in the platform local reference frame after
            applying the platform_bearing rotation.
        """
        rotated_offsets = []
        for offset in mooring_offsets:
            x = offset[0]*np.cos(np.deg2rad(self.platform_bearing)) + offset[1]*np.sin(np.deg2rad(self.platform_bearing))
            y = -offset[0]*np.sin(np.deg2rad(self.platform_bearing)) + offset[1]*np.cos(np.deg2rad(self.platform_bearing))
            rotated_offsets.append((x, y))
        return rotated_offsets
    
    def _set_mooring_bearings(self, mooring_bearings: tuple[float, ...] | None) -> None:
        """
        Private method to set the mooring Cardinal Bearings after the platform bearing rotation.

        Parameters
        ----------
        mooring_bearings : tuple[float, ...] | None
            Mooring cardinal bearings in the platform local frame of reference before applying the
            platform_bearing rotation.

        Raises
        ------
        IndexError
            The number of moorings does not match the number of mooring bearings pairs passed.
        """
        if mooring_bearings is None:
            if self.n_moorings % 2:
                self.mooring_bearings = (np.arange(0, -360, -360/self.n_moorings) + self.platform_bearing) % 360
            else:
                self.mooring_bearings = (np.arange(360/self.n_moorings/2, -(360 - 360/self.n_moorings/2), -360/self.n_moorings) + self.platform_bearing) % 360
        elif len(mooring_bearings) != self.n_moorings:
            raise IndexError(f"length of mooring_bearings must match n_moorings. {mooring_bearings=} and {self.n_moorings=} passed.")
        else:
            self.mooring_bearings = (np.array(mooring_bearings) + self.platform_bearing) % 360
    
    def _calculate_mooring_origins(self) -> None:
        """
        Private method to calculate all mooring fairlead coordinates i nthe global frame of \
        reference for each turbine platform.
        """
        mooring_origins = []
        for point in self.geom.geoms:
            local_origins = []
            for offset in self.mooring_offsets:
                local_origins.append((offset[0] + point.x, offset[1] + point.y))
            mooring_origins.append(local_origins)
        
        self.mooring_origins = mooring_origins
    
    def _calculate_mooring_radii(self) -> None:
        """
        Private method to calculate the mooring radius at each turbine location.

        Raises
        ------
        AttributeError
            If the mooring radius function is a Callable and no bathymetry dataset is loaded for
            the layout.
        TypeError
            If the mooring radius function is not a float or a Callable.
        """
        if callable(self.mooring_radius_func):
            if self.samples is None:
                raise AttributeError("No bathymetry loaded for the layout.")
            else:
                bathymetry = np.array(self.samples)
                bathy_func = np.vectorize(self.mooring_radius_func)
                radii = bathy_func(bathymetry)
                self.mooring_radii = radii
        elif isinstance(self.mooring_radius_func, float) or isinstance(self.mooring_radius_func, int):
            radii = np.full(self.n_wtg, fill_value=self.mooring_radius_func)
            self.mooring_radii = radii
        else:
            raise TypeError(f"Invalid type for radius. Must be float or callable. Type of {type(self.mooring_radius_func)} passed.")
    
    def _calc_mooring_anchors(self) -> None:
        mooring_anchors = []
        for platform_mooring_origins, radius in zip(self.mooring_origins, self.mooring_radii):
            platform_anchors = []
            for i, origin in enumerate(platform_mooring_origins):
                x = origin[0] + radius*np.sin(np.deg2rad(self.mooring_bearings[i]))
                y = origin[1] + radius*np.cos(np.deg2rad(self.mooring_bearings[i]))
                anchor = (x, y)
                platform_anchors.append(anchor)
            mooring_anchors.append(platform_anchors)
        
        self.mooring_anchors = mooring_anchors
    
    def _generate_mooring_geometry(self) -> None:
        """
        Private method to generate the mooring MultiLineString geometries.
        """
        self._calculate_mooring_origins()
        self._calculate_mooring_radii()
        self._calc_mooring_anchors()
        
        mooring_geom = []
        
        for origins, anchors in zip(self.mooring_origins, self.mooring_anchors):
            mooring_geom.append(MultiLineString(list(zip(origins, anchors))))
        
        self.mooring_geom = mooring_geom
    
    def generate_moorings(self,
                          n_moorings: int,
                          radius: float | Callable[[float], float],
                          platform_bearing: float | None = None,
                          mooring_bearings: tuple[float, ...] | None = None,
                          mooring_offsets: list[tuple[float, float]] | None = None,
                          contain: bool = False,
                          mooring_areas: list[Polygon | MultiPolygon] | None = None,
                          avoid_exclusions: bool = False,
                          mooring_exclusions: list[Polygon | MultiPolygon] | None = None
                          ) -> None:
        """
        This method is used to generate mooring MultiLineString geometries for the FloatingLayout \
        instance.

        Parameters
        ----------
        n_moorings : int
            Number of mooring lines at each floating platfrom.
        radius : float | Callable[[float], float]
            Mooring radius function. This parameter defines how the mooring radius length is
            calculated for each turbine floating platform. If a float is passed, it is used as a
            constant mooring length for all turbine positions. If a callable is passed, it must
            take in a single argument of bathymetry for each turbine postion.
        platform_bearing : float | None, optional
            Cardinal bearing of the floating platform with respect to True North, by default None.
        mooring_bearings : tuple[float, ...] | None, optional
            Cardinal bearings of the mooring lines in the local reference frame of the floating
            platform. If None is passed, the moorings will be spread evenly around the 360°
            anti-clockwise, by default None.
        mooring_offsets : list[tuple[float, float]] | None, optional
            The coordinate pairs for the mooring fairleads with respect to the turbine position in
            the floating platform local reference frame. If None is passed, all mooring fairleads
            will be at the turbine location, by default None.
        contain : bool, optional
            If True, the mooring lines and associated turbine locations outwith the mooring_areas
            are removed from the layout. If set to True but no mooring_areas are passed, the
            normal turbine area constraint geometry will be used, by default False.
        mooring_areas : list[Polygon  |  MultiPolygon] | None, optional
            Polygon geometries to constrain the mooring lines to, used in conjunction with the
            contain argument, by default None. Has no effect if contain is set to False. 
        avoid_exclusions : bool, optional
            If True, removes the mooring lines and associated turbine locations that overlap with
            the mooring_exclusions. If set to True but no mooring_exclusions are passed, the normal
            turbine exclusion geometry will be used, by default False.
        mooring_exclusions : list[Polygon  |  MultiPolygon] | None, optional
            Polygon geometries to use as mooring exclusion zones, used in conjunction with the
            avoid_exclusions argument, by default None. Has no effect if avoid_exclusions is set to False.

        Raises
        ------
        AttributeError
            If the FloatingLayout instance has no turbine geometries from which to generate
            moorings from.
        """
        if not self.has_geom():
            raise AttributeError("Layout has no point geometry to generate moorings from.")

        self.n_moorings = n_moorings
        self.mooring_radius_func = radius
        self._set_platform_bearing(platform_bearing)
        self._set_mooring_offsets(mooring_offsets)
        self._set_mooring_bearings(mooring_bearings)
            
        self._generate_mooring_geometry()
        
        if contain:
            self.clip_moorings_to_area(areas=mooring_areas, mode='w')
        if avoid_exclusions:
            self.avoid_mooring_exclusions(exclusions=mooring_exclusions, mode='w')
    
    def reset_moorings(self) -> "FloatingLayout":
        """
        This method resets all moorings and any mooring related attributes back to None. \
        Also recovers the turbine positions lost due to any mooring contraints.

        Returns
        -------
        FloatingLayout
            FloatingLayout with reset moorings and associated attributes.
        """
        self.mooring_geom = None
        self.n_moorings = None
        self.mooring_radius_func = None
        self.platform_bearing = None
        self.mooring_bearings = None
        self.mooring_offsets = None
        self.mooring_origins = None
        self.mooring_anchors = None
        self.mooring_radii = None
        self.mooring_area = None
        self.mooring_exclusions = None
        self.mooring_buffers = None
        
        if self.has_geom():
            self.geom = self._raw_geom
            if self.has_area():
                self.clip_to_area(self.area)
            if self.has_exclusion():
                self.avoid_exclusions(self.exclusion)
            if self.has_bathymetry():
                self.load_bathymetry(
                                    filepath=self.bathymetry_path,
                                    limits=self.bathymetry_limits,
                                    sign=self.bathymetry_sign,
                                    drop_na=self.bathymetry_drop_na)
        return self
    
    @property
    def has_moorings(self) -> bool:
        """
        Returns True if the FloatingLayout has mooring geometries.

        Returns
        -------
        bool
            True if the FloatingLayout has mooring geometries, False otherwise.
        """
        return self.mooring_geom is not None
    
    @property
    def has_mooring_area(self) -> bool:
        """
        Returns True if the FloatingLayout has mooring areas.

        Returns
        -------
        bool
            True if the FloatingLayout has mooring areas, False otherwise.
        """
        return self.mooring_area is not None
    
    @property
    def has_mooring_exclusions(self) -> bool:
        """
        Returns True if the FloatingLayout has mooring exclusions.

        Returns
        -------
        bool
            True if the FloatingLayout has mooring exclusions, False otherwise.
        """
        return self.mooring_exclusions is not None
    
    def reset_mooring_area(self) -> "FloatingLayout":
        """
        Resets the mooring_area constraints and recovers any turbine positions lost because of it.

        Returns
        -------
        FloatingLayout
            FloatingLayout without any mooring_area constraints.
        """
        if self.has_mooring_area:
            self.mooring_area = None
            if self.has_geom():
                self.geom = self._raw_geom
                if self.has_area():
                    self.clip_to_area(self.area)
                if self.has_exclusion():
                    self.avoid_exclusions(self.exclusion)
                if self.has_bathymetry():
                    self.load_bathymetry(
                                        filepath=self.bathymetry_path,
                                        limits=self.bathymetry_limits,
                                        sign=self.bathymetry_sign,
                                        drop_na=self.bathymetry_drop_na)
                
                self.generate_moorings(n_moorings=self.n_moorings,
                                       radius=self.mooring_radius_func,
                                       platform_bearing=self.platform_bearing,
                                       mooring_bearings=self.mooring_bearings,
                                       mooring_offsets=self.mooring_offsets
                                       )
                if self.mooring_exclusions is not None:
                    self.avoid_mooring_exclusions(exclusions=[self.mooring_exclusions], mode='w')
        return self
    
    def reset_mooring_exclusions(self) -> "FloatingLayout":
        """
        Resets the mooring_exclusions constraints and recovers any turbine positions lost because of it.

        Returns
        -------
        FloatingLayout
            FloatingLayout without any mooring_exclusions constraints.
        """
        if self.has_mooring_exclusions:
            self.mooring_exclusions = None
            if self.has_geom():
                self.geom = self._raw_geom
                if self.has_area():
                    self.clip_to_area(self.area)
                if self.has_exclusion():
                    self.avoid_exclusions(self.exclusion)
                if self.has_bathymetry():
                    self.load_bathymetry(
                                        filepath=self.bathymetry_path,
                                        limits=self.bathymetry_limits,
                                        sign=self.bathymetry_sign,
                                        drop_na=self.bathymetry_drop_na)
                
                self.generate_moorings(n_moorings=self.n_moorings,
                                       radius=self.mooring_radius_func,
                                       platform_bearing=self.platform_bearing,
                                       mooring_bearings=self.mooring_bearings,
                                       mooring_offsets=self.mooring_offsets
                                       )
                if self.mooring_area is not None:
                    self.clip_moorings_to_area(areas=[self.mooring_area], mode='w')
        return self
            
    
    def clip_moorings_to_area(self, areas: list[Polygon | MultiPolygon] | None = None, mode: str = "a") -> "FloatingLayout":
        """
        Applies a mooring_area constraint to the floating Layout. Moorings and associated turbine \
        positions that fall outside the area are removed.

        Parameters
        ----------
        areas : list[Polygon  |  MultiPolygon] | None, optional
            Polygoin geometries to use as mooring_area constraints. If None is passed, and the
            FloatingLayout has an area constraint, the normal area constraint will be used, by default None.
        mode : str, optional
            Speficies the mode in which the passed areas are added to the layout as
            constraints. Options are 'w' for write mode and 'a' for append mode (case insensitive).
            When calling this method in append mode, the passed geometries will be added to any
            existing mooring_areas as additional constraints. When calling this method in write
            mode, any existing area constraints will be overwritten by the passed geometries. By
            default "a".

        Returns
        -------
        FloatingLayout
            FloatingLayout with mooring lines, and associated turbines, limited to within the mooring_areas constraint.

        Raises
        ------
        AttributeError
            If the FloatingLayout has no mooring geometries to clip.
        ValueError
            If mode is not 'a' or 'w'.
        AttributeError
            If no areas are passed and the FloatingLayout does not have an area constraint.
        """
        if not self.has_moorings:
            raise AttributeError(f"{self.__class__} has no moorings to clip.")
        
        if mode.lower() not in ["a", "w"]:
            raise ValueError(f"mode must be 'a' or 'w'. Value of {mode=} passed.")
        
        if areas is None:
            if not self.has_area():
                raise AttributeError(f"No area geometries passed and {self.__class__} has no areas defined.")
            else:
                self.mooring_area = self.area
        else:
            area = unary_union(areas)
            if self.has_mooring_area and mode.lower() == 'a':
                self.mooring_area = unary_union([self.mooring_area, area])
            elif self.has_mooring_area and mode.lower() == 'w':
                self.reset_mooring_area()
                self.mooring_area = area
            else:
                self.mooring_area = area
        
        valid_points = []
        valid_moorings = []
        valid_mooring_origins = []
        valid_mooring_anchors = []
        valid_mooring_radii = []
        if self.has_bathymetry():
            valid_samples = []
        
        for i, mooring in enumerate(self.mooring_geom):
            if mooring.within(self.mooring_area):
                valid_points.append(self.geom.geoms[i])
                valid_mooring_origins.append(self.mooring_origins[i])
                valid_mooring_anchors.append(self.mooring_anchors[i])
                valid_mooring_radii.append(self.mooring_radii[i])
                valid_moorings.append(mooring)
                if self.has_bathymetry():
                    valid_samples.append(self.samples[i])
        
        self.geom = unary_union(valid_points)
        self.mooring_origins = valid_mooring_origins
        self.mooring_anchors = valid_mooring_anchors
        self.mooring_radii = np.array(valid_mooring_radii)
        self.mooring_geom = valid_moorings
        if self.has_bathymetry():
            self.samples = valid_samples
        
        return self

    def avoid_mooring_exclusions(self, exclusions: list[Polygon | MultiPolygon] | None = None, mode: str = "a") -> "FloatingLayout":
        """
        Applies a mooring_exclusions constraint to the floating Layout. Moorings and associated \
        turbine positions that overlap with these exclusions are removed fro mthe layout.

        Parameters
        ----------
        exclusions : list[Polygon  |  MultiPolygon] | None, optional
            Polygoin geometries to use as mooring_exclusion constraints. If None is passed, and the
            FloatingLayout has an exclusion constraint, the normal exclusion constraint will be
            used, by default None.
        mode : str, optional
            Speficies the mode in which the passed exclusions are added to the layout as
            constraints. Options are 'w' for write mode and 'a' for append mode (case insensitive).
            When calling this method in append mode, the passed geometries will be added to any
            existing moring_exclusions as additional constraints. When calling this method in write
            mode, any existing mooring_exclusion constraints will be overwritten by the passed geometries. By
            default "a".

        Returns
        -------
        FloatingLayout
            FloatingLayout with no mooring lines and associated turbines overlapping any
            mooring_exclusions.

        Raises
        ------
        AttributeError
            If the FloatingLayout has no mooring geometries.
        ValueError
            If mode is not 'a' or 'w'.
        AttributeError
            If no exclusions are passed and the FloatingLayout does not have an exclusion constraint.
        """
        if not self.has_moorings:
            raise AttributeError(f"{self.__class__} has no moorings.")
        
        if mode.lower() not in ["a", "w"]:
            raise ValueError(f"mode must be 'a' or 'w'. Value of {mode=} passed.")
        
        if exclusions is None:
            if not self.has_exclusion():
                raise AttributeError(f"No exclusions geometries passed and {self.__class__} has no exclusions defined.")
            else:
                self.mooring_exclusions = self.exclusion
        else:
            exclusion = unary_union(exclusions)
            if self.has_mooring_exclusions and mode.lower() == 'a':
                self.mooring_exclusions = unary_union([self.mooring_exclusions, exclusion])
            elif self.has_mooring_exclusions and mode.lower() == 'w':
                self.reset_mooring_exclusions()
                self.mooring_exclusions = exclusion
            else:
                self.mooring_exclusions = exclusion
        
        valid_points = []
        valid_moorings = []
        valid_mooring_origins = []
        valid_mooring_anchors = []
        valid_mooring_radii = []
        if self.has_bathymetry():
            valid_samples = []
        
        for i, mooring in enumerate(self.mooring_geom):
            if mooring.disjoint(self.mooring_exclusions):
                valid_points.append(self.geom.geoms[i])
                valid_mooring_origins.append(self.mooring_origins[i])
                valid_mooring_anchors.append(self.mooring_anchors[i])
                valid_mooring_radii.append(self.mooring_radii[i])
                valid_moorings.append(mooring)
                if self.has_bathymetry():
                    valid_samples.append(self.samples[i])
        
        self.geom = unary_union(valid_points)
        self.mooring_origins = valid_mooring_origins
        self.mooring_anchors = valid_mooring_anchors
        self.mooring_radii = np.array(valid_mooring_radii)
        self.mooring_geom = valid_moorings
        if self.has_bathymetry():
            self.samples = valid_samples
        
        return self
    
    @property
    def mooring_crossings(self) -> bool:
        """
        Returns True if there are any mooring crossings in the x-y plane.

        Returns
        -------
        bool
            True if any mooring line geometries cross in the x-y plane. False otherwise.
        """
        if not self.has_moorings:
            return False
        
        for i, mooring in enumerate(self.mooring_geom):
            if i == 0:
                other_moorings = unary_union(self.mooring_geom[1:])
            elif i == len(self.mooring_geom) - 1:
                other_moorings = unary_union(self.mooring_geom[:-1])
            else:
                other_moorings = unary_union(self.mooring_geom[0:i] + self.mooring_geom[i+1:])
    
            if mooring.intersects(other_moorings):
                return True
        
        return False
    
    def generate_mooring_buffers(self, width: float = 10, segments: int =5) -> None:
        """
        Generates angular buffers around the layout mooring lines.

        Parameters
        ----------
        width : float, optional
            Angular width of the buffer around the mooring lines in degrees. Creates a buffer of
            half the width to either side of the mooring lines. By default 10.
        segments : int, optional
            Number of segments to use to approximate the circular arc, by default 5.

        Raises
        ------
        AttributeError
            If the layout has no mooring geometries to generate buffers around.
        """
        if not self.has_moorings:
            raise AttributeError("Layout does not have any moorings.")

        buffers = []

        for xy_wtg, r in zip(self.coords, self.mooring_radii):
            for bearing, offset in zip(self.mooring_bearings, np.array(self.mooring_offsets)):
                start_angle = bearing - width / 2
                end_angle = bearing + width / 2
                
                xy = np.empty((segments + 2, 2))
                xy[0] = xy_wtg + offset
                theta = np.radians(np.linspace(start_angle, end_angle, segments + 1))
                xy[1:] = xy_wtg + offset + r * np.stack([np.sin(theta), np.cos(theta)], axis=-1)
                buffers.append(Polygon(xy))
        
        self.mooring_buffers = MultiPolygon(buffers)
    
    @property
    def has_mooring_buffers(self) -> bool:
        """
        Returns True if there are any generate mooring buffers.

        Returns
        -------
        bool
            True if mooring buffers have been generated. False otherwise.
        """
        return self.mooring_buffers is not None

    @property
    def mooring_buffer_overlap(self) -> bool:
        """
        Returns True if the mooring_buffers overlap.

        Returns
        -------
        bool
            True if mooring_buffers overlap, False otherwise.
        """
        if not self.has_mooring_buffers:
            return False
        
        n_wtg = self.n_wtg
        n_moorings = self.n_moorings
        
        for i in range(n_wtg):
            buffer = self.mooring_buffers.geoms[i*n_moorings:i*n_moorings + n_moorings]
            if i == 0:
                other_buffers = unary_union(self.mooring_buffers.geoms[i*n_moorings + n_moorings:])
            elif i == n_wtg - 1:
                other_buffers = unary_union(self.mooring_buffers.geoms[:-n_moorings])
            else:
                other_buffers = unary_union([self.mooring_buffers.geoms[0:i*n_moorings], self.mooring_buffers.geoms[i*n_moorings + n_moorings:]])
    
            if buffer.intersects(other_buffers):
                return True
        
        return False

    def plot(self, ax: Axes | None = None, figsize: tuple = (8, 6), show_bathy: bool = False) -> tuple[Figure, Axes]:
        """
        Plots the FloatingLayout turbine and mooring geometries, as well as areas and exclusions.

        Parameters
        ----------
        ax : matplotlib.axes.Axes | None, optional
            Axes on which the Layout is drawn. If None, a new Figure and Axes are created, by default None.
        figsize : tuple, optional
            Figure size in inches in (width, height) format, by default (8, 6)
        show_bathy : bool, optional
            If True, render the bathymetry onto the plot. Layout must have a loaded bathymetry filepath.

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
            Tuple of Figure and Axes objects.

        Raises
        ------
        TypeError
            The passed ax argument was not of type Axes.
        
        Caution
        ----
        Performance can be seriously reduced if show_bathy is set to true and high resolution bathymetry data is to be rendered.
        """
        fig, ax = super().plot(ax=ax, figsize=figsize, show_bathy=show_bathy)
        if self.mooring_geom is not None:
            gp.GeoSeries(self.mooring_geom).plot(ax=ax, color='black')
        
        if self.has_mooring_area:
            gp.GeoSeries(self.mooring_area).plot(ax=ax, facecolor='green', alpha=0.2)
        
        if self.has_mooring_exclusions:
             gp.GeoSeries(self.mooring_exclusions).plot(ax=ax, facecolor='red', alpha=0.2)
        
        if self.has_mooring_buffers:
            gp.GeoSeries(self.mooring_buffers).plot(ax=ax, facecolor='grey', edgecolor='black', hatch="//", alpha=0.4)
        
        return fig, ax
    
    def to_shapefile(self, filepath: str, epsg: int = None) -> None:
        """
        Saves the FloatingLayout turbine positions as Point geometries in a Shapefile.\
        If the FloatingLayout has moorings and mooring buffers, the moorings are saved as MultiLine \
        geometries and the mooring buffers as Polygon in separate Shapefiles.

        Parameters
        ----------
        filepath : str
            Filepath to save the Shapefile as. The mooring shapefile will have the "moorings_" prefix.
        epsg : int
            Valid `EPSG <https://epsg.io/>`_ code for the projection of the Shapefile(s). If passed,
            this code takes precedence over the layout epsg code (if defined). By default, None.
        
        Raises
        ------
        ValueError
            No EPSG code available for the shapefile CRS.
        """
        
        super().to_shapefile(filepath, epsg)
        
        if epsg is None and self.epsg is not None:
            crs = pyproj.CRS.from_epsg(self.epsg).to_dict()
        elif epsg is None and self.epsg is None:
            raise ValueError("No EPSG code provided.")
        elif _is_epsg(epsg):
            crs = pyproj.CRS.from_epsg(epsg).to_dict()
        
        if self.has_moorings:
            # TODO: Create a more detailed attribute schema (wakes, water depth, etc.).
            schema = {"geometry": "MultiLineString", "properties": {"id": "int", "radius": "float"}}

            driver = "ESRI Shapefile"
            
            if '.' in filepath:
                mooring_path = "_moorings.".join(filepath.split('.'))
            else:
                mooring_path = filepath + "_moorings"

            with fiona.open(
                mooring_path, "w", crs=crs, driver=driver, schema=schema
            ) as new_file:
                for idx, geom in enumerate(self.mooring_geom, start=1):  # type: ignore
                    new_file.write({"geometry": mapping(geom), "properties": {"id": idx, "radius": float(self.mooring_radii[idx - 1])}})
            
            if self.has_mooring_buffers:
                schema = {"geometry": "Polygon", "properties": {"id": "int"}}
                driver = "ESRI Shapefile"
                
                if '.' in filepath:
                    cone_path = "_mooring_buffers.".join(filepath.split('.'))
                else:
                    cone_path = filepath + "_mooring_buffers"
                
                with fiona.open(
                    cone_path, "w", crs=crs, driver=driver, schema=schema
                ) as new_file:
                    for idx, geom in enumerate(self.mooring_buffers.geoms, start=1):  # type: ignore
                        new_file.write({"geometry": mapping(geom), "properties": {"id": idx}})

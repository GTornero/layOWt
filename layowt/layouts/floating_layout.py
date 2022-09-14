"""
This module defines the FloatingLayout class.
This class can be used to design and represent floating turbines and their mooring lines.
"""

from typing import Callable

import numpy as np
from layowt.grids import Grid
from layowt.layouts.layout import Layout
from shapely.geometry import MultiLineString, MultiPolygon, Polygon


class FloatingLayout(Layout):
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
    
    def _set_platform_bearing(self, platform_bearing: float | None) -> None:
        if platform_bearing is None:
            self.platform_bearing = 0
        else:
            self.platform_bearing = platform_bearing
    
    def _set_mooring_offsets(self, mooring_offsets: list[tuple[float, float]] | None) -> None:
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
        rotated_offsets = []
        for offset in mooring_offsets:
            x = offset[0]*np.cos(np.deg2rad(self.platform_bearing)) + offset[1]*np.sin(np.deg2rad(self.platform_bearing))
            y = -offset[0]*np.sin(np.deg2rad(self.platform_bearing)) + offset[1]*np.cos(np.deg2rad(self.platform_bearing))
            rotated_offsets.append((x, y))
        return rotated_offsets
    
    def _set_mooring_bearings(self, mooring_bearings: tuple[float, ...] | None) -> None:
        if mooring_bearings is None:
            if self.n_moorings % 2:
                self.mooring_bearings = (np.arange(0, -360, -360/self.n_moorings) + self.platform_bearing) % 360
            else:
                self.mooring_bearings = (np.arange(360/self.n_moorings/2, -(360 - 360/self.n_moorings/2), -360/self.n_moorings) + self.platform_bearing) % 360
        elif len(mooring_bearings) != self.n_moorings:
            raise IndexError(f"length of mooring_bearings must match n_moorings. {mooring_bearings=} and {self.n_moorings=} passed.")
        else:
            self.mooring_bearings = (np.array(mooring_bearings) + self.platform_bearing) % 360
    
    def _calculate_mooring_origins(self) -> list[list[tuple[float, float]]]:
        mooring_origins = []
        for point in self.geom.geoms:
            local_origins = []
            for offset in self.mooring_offsets:
                local_origins.append((offset[0] + point.x, offset[1] + point.y))
            mooring_origins.append(local_origins)
        
        return mooring_origins
    
    def _calculate_mooring_radii(self) -> np.ndarray:
        if callable(self.mooring_radius_func):
            if self.samples is None:
                raise AttributeError("No bathymetry loaded for the layout.")
            else:
                bathymetry = np.array(self.samples)
                bathy_func = np.vectorize(self.mooring_radius_func)
                radii = bathy_func(bathymetry)
                return radii
        elif isinstance(self.mooring_radius_func, float) or isinstance(self.mooring_radius_func, int):
            radii = np.full(self.n_wtg, fill_value=self.mooring_radius_func)
            return radii
        else:
            raise TypeError(f"Invalid type for radius. Must be float or callable. Type of {type(self.mooring_radius_func)} passed.")
    
    def _calc_mooring_anchors(self) -> list[list[tuple[float, float]]]:
        mooring_anchors = []
        for platform_mooring_origins, radius in zip(self.mooring_origins, self.mooring_radii):
            platform_anchors = []
            for i, origin in enumerate(platform_mooring_origins):
                x = origin[0] + radius*np.sin(np.deg2rad(self.mooring_bearings[i]))
                y = origin[1] + radius*np.cos(np.deg2rad(self.mooring_bearings[i]))
                anchor = (x, y)
                platform_anchors.append(anchor)
            mooring_anchors.append(platform_anchors)
        
        return mooring_anchors
    
    def _generate_mooring_geometry(self) -> None:
        
        self.mooring_origins = self._calculate_mooring_origins()
        self.mooring_radii = self._calculate_mooring_radii()
        self.mooring_anchors = self._calc_mooring_anchors()
        
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
                          allow_crossings: bool = True,
                          contain: bool = False,
                          avoid_exclusions: bool = False,
                          exclusions: list[Polygon | MultiPolygon] | None = None
                          ) -> None:
        # TODO: refactor this into a MooringLayout class within the Floating Layout class.
        self.n_moorings = n_moorings
        self.mooring_radius_func = radius
        # TODO:
            # 1: set platform bearing
        
        self._set_platform_bearing(platform_bearing)
            # 2: set local mooring offsets
        self._set_mooring_offsets(mooring_offsets)
            # 3: calculate mooring bearings / set mooring bearings
        self._set_mooring_bearings(mooring_bearings)
            # for each wtg position:
        self._generate_mooring_geometry()
                # 4: calculate mooring origin location
                # 5: calculate mooring radius
                # 6: calculate anchor location
                # 7: generate mooring geometry (list of lines or multiline??)
                # 8: (optional) eliminate positions with moorings outside lease area
                # 9: (optional) eliminate positions crossing exclusions (can be another exclusion layer)
                # 10: (optional) eliminate positions with crossing moorings (which of the two positions involved in the crossing?)
            

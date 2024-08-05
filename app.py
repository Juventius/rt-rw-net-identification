import streamlit as st
import pandas as pd
from scipy.spatial import KDTree
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from fastkml import kml
from io import BytesIO

# Define the earth radius in kilometers
EARTH_RADIUS_KM = 6371.0

# Function to convert degrees to radians
def deg2rad(deg):
    return deg * (np.pi / 180)

# Function to calculate the distance between two points (Haversine formula)
def haversine(lat1, lon1, lat2, lon2):
    dlat = deg2rad(lat2 - lat1)
    dlon = deg2rad(lat2 - lon2)
    a = np.sin(dlat / 2) ** 2 + np.cos(deg2rad(lat1)) * np.cos(deg2rad(lat2)) * np.sin(dlon / 2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    return EARTH_RADIUS_KM * c

st.title("RT/RW Net Identification")
st.markdown("made by Juventius Kriswijanarko \([LinkedIn](https://www.linkedin.com/in/juventius-kriswijanarko/)\)")
st.markdown("*Please contact the above person for feedbacks and suggestions.")
st.markdown("*[Sample Datasets](https://drive.google.com/drive/folders/1loDyXbKmJ31iNdCnSLI1wIkhD_8nstZS?usp=sharing) (current access for erickhimura@gmail.com)")

ookla_file = st.file_uploader("Upload Ookla Excel Dataset:", type=["xlsx"])
odp_file = st.file_uploader("Upload ODP Excel Dataset:", type=["xlsx"])

# Add a slider for the user to set the latency threshold
latency_threshold = st.slider("Set the latency threshold (in ms):", min_value=0, max_value=200, value=55)

if 'processed' not in st.session_state:
    st.session_state.processed = False

if st.button("Process Files"):
    if ookla_file and odp_file:
        st.write("Processing, please wait... (this might take approx. 5-10 minutes)")

        # Read the Ookla file
        st.write("Reading Ookla file...")
        signal_df = pd.read_excel(ookla_file, engine="openpyxl")
        st.write("Done.")

        # Read the ODP file
        st.write("Reading ODP file...")
        odp_df = pd.read_excel(odp_file, engine="openpyxl")
        st.write("Done.")

        # Filtering and Cleaning Ookla Dataset
        st.write("Filtering and cleaning Ookla dataset...")
        st.write('Checking attributes: "DESA", "KECAMATAN", "KABUPATEN", "provider", "attr_location_latitude", "attr_location_longitude", "val_latency_iqm_ms"')
        filtered_signal_df = signal_df[["DESA", "KECAMATAN", "KABUPATEN", "provider", "attr_location_latitude", "attr_location_longitude", "val_latency_iqm_ms"]]
        cleaned_signal_df = filtered_signal_df.dropna(subset=["val_latency_iqm_ms"])
        cleaned_signal_df.reset_index(drop=True, inplace=True)
        st.write("Clear.")

        # Filtering ODP Dataset
        st.write("Filtering ODP dataset...")
        st.write('Checking attributes: "ODP NAME", "LATITUDE", "LONGITUDE"')
        filtered_odp_df = odp_df[["ODP NAME", "LATITUDE", "LONGITUDE"]]
        filtered_odp_df.reset_index(drop=True, inplace=True)
        st.write("Clear.")

        # Filter the signal points with val_latency_iqm_ms over latency_threshold
        filtered_signals = cleaned_signal_df[cleaned_signal_df["val_latency_iqm_ms"] > latency_threshold]

        # Create KDTree for ODP locations
        odp_locations = filtered_odp_df[["LATITUDE", "LONGITUDE"]].values
        odp_tree = KDTree(odp_locations)

        # Initialize the result lists
        signal_points_for_improvement = []
        potential_rt_rw_net = []

        # Analysis
        st.write("Analyzing...")
        for _, signal in filtered_signals.iterrows():
            signal_desa = signal["DESA"]
            signal_kecamatan = signal["KECAMATAN"]
            signal_kabupaten = signal["KABUPATEN"]
            signal_provider = signal["provider"]
            signal_lat = signal["attr_location_latitude"]
            signal_lon = signal["attr_location_longitude"]
            signal_latency = signal["val_latency_iqm_ms"]
            
            # Search for ODPs within 500 meters and 1 kilometer
            odp_in_500m = False
            odp_in_1km = False

            # Radius in km
            radius_500m = 0.5
            radius_1km = 1.0

            indices_500m = odp_tree.query_ball_point([signal_lat, signal_lon], radius_500m / EARTH_RADIUS_KM * (180 / np.pi))
            indices_1km = odp_tree.query_ball_point([signal_lat, signal_lon], radius_1km / EARTH_RADIUS_KM * (180 / np.pi))

            if indices_500m:
                odp_in_500m = True
            
            if indices_1km:
                odp_in_1km = True

            if not odp_in_500m and odp_in_1km:
                signal_points_for_improvement.append({
                    "desa": signal_desa,
                    "kecamatan": signal_kecamatan,
                    "kabupaten": signal_kabupaten,
                    "provider": signal_provider,
                    "latitude": signal_lat,
                    "longitude": signal_lon,
                    "latency": signal_latency
                })
            elif not odp_in_500m and not odp_in_1km:
                potential_rt_rw_net.append({
                    "desa": signal_desa,
                    "kecamatan": signal_kecamatan,
                    "kabupaten": signal_kabupaten,
                    "provider": signal_provider,
                    "latitude": signal_lat,
                    "longitude": signal_lon,
                    "latency": signal_latency
                })

        telkomsel_providers = ["PT Telekomunikasi Indonesia", "Telkomsel Orbit", "by.U"]

        # Generate dataframes
        signal_points_for_improvement_df = pd.DataFrame(signal_points_for_improvement)
        signal_points_for_improvement_df = signal_points_for_improvement_df[signal_points_for_improvement_df["provider"].isin(telkomsel_providers)]
        signal_points_for_improvement_df.reset_index(drop=True, inplace=True)
        potential_rt_rw_net_df = pd.DataFrame(potential_rt_rw_net)

        # Split RT/RW Net Dataframe
        telkomsel_rt_rw_net_df = potential_rt_rw_net_df[potential_rt_rw_net_df["provider"].isin(telkomsel_providers)]
        telkomsel_rt_rw_net_df.reset_index(drop=True, inplace=True)
        others_rt_rw_net_df = potential_rt_rw_net_df[~potential_rt_rw_net_df["provider"].isin(telkomsel_providers)]
        others_rt_rw_net_df.reset_index(drop=True, inplace=True)

        st.write("All processes are finished.")
        st.write("===========================================")
        st.session_state.signal_points_for_improvement_df = signal_points_for_improvement_df
        st.session_state.telkomsel_rt_rw_net_df = telkomsel_rt_rw_net_df
        st.session_state.others_rt_rw_net_df = others_rt_rw_net_df
        st.session_state.filtered_odp_df = filtered_odp_df
        st.session_state.processed = True

if st.session_state.processed:
    st.write("1. Signal points for improvement (500 meters):")
    st.dataframe(st.session_state.signal_points_for_improvement_df)
    st.write("2. Potential RT/RW Net points (Telkomsel):")
    st.dataframe(st.session_state.telkomsel_rt_rw_net_df)
    st.write("3. Potential RT/RW Net points (Others):")
    st.dataframe(st.session_state.others_rt_rw_net_df)

    export_option = st.selectbox("Which file do you want to export?", ["Select", "ODP", "500m", "RT/RW (Telkomsel)", "RT/RW (Others)"])

    if export_option == "ODP":
        # Convert filtered_odp_df into a GeoDataFrame
        st.write("Exporting ODP file...")
        geometry = [Point(xy) for xy in zip(st.session_state.filtered_odp_df["LONGITUDE"], st.session_state.filtered_odp_df["LATITUDE"])]
        geo_df = gpd.GeoDataFrame(st.session_state.filtered_odp_df, geometry=geometry, crs="EPSG:4326")
        
        # Save the GeoDataFrame to a GeoPackage
        geopackage_buffer = BytesIO()
        geo_df.to_file(geopackage_buffer, driver="GPKG", layer='ODP')
        geopackage_buffer.seek(0)
        st.write("Done.")

        st.download_button(
            label="Download ODP GeoPackage",
            data=geopackage_buffer,
            file_name="ODP.gpkg",
            mime="application/octet-stream"
        )
    
    elif export_option == "500m":
        # Convert signal_points_for_improvement_df into a GeoDataFrame
        st.write("Exporting Radius 500m file...")
        geometry = [Point(xy) for xy in zip(st.session_state.signal_points_for_improvement_df["longitude"], st.session_state.signal_points_for_improvement_df["latitude"])]
        geo_df = gpd.GeoDataFrame(st.session_state.signal_points_for_improvement_df, geometry=geometry, crs="EPSG:4326")
        
        # Save the GeoDataFrame to a GeoPackage
        geopackage_buffer = BytesIO()
        geo_df.to_file(geopackage_buffer, driver="GPKG", layer='Radius_500m')
        geopackage_buffer.seek(0)
        st.write("Done.")

        st.download_button(
            label="Download Radius 500m GeoPackage",
            data=geopackage_buffer,
            file_name="Radius_500m.gpkg",
            mime="application/octet-stream"
        )
    
    elif export_option == "RT/RW (Telkomsel)":
        # Initialize KML object
        st.write("Exporting RT/RW Net (Telkomsel) file...")
        k = kml.KML()
        
        # Create a KML document
        document = kml.Document()
        k.append(document)
        
        # Create a KML placemark for each point
        for _, row in st.session_state.telkomsel_rt_rw_net_df.iterrows():
            placemark = kml.Placemark()
            placemark.geometry = Point(row["longitude"], row["latitude"])
            placemark.name = f"{row['provider']} (Lat: {row['latitude']}, Lon: {row['longitude']})"
            
            # Add attributes as KML extended data
            extended_data = kml.ExtendedData()
            extended_data.elements.append(kml.Data(name="desa", value=str(row["desa"])))
            extended_data.elements.append(kml.Data(name="kecamatan", value=str(row["kecamatan"])))
            extended_data.elements.append(kml.Data(name="kabupaten", value=str(row["kabupaten"])))
            extended_data.elements.append(kml.Data(name="provider", value=str(row["provider"])))
            extended_data.elements.append(kml.Data(name="latency", value=str(row["latency"])))
            placemark.extended_data = extended_data
            
            document.append(placemark)
        
        # Save the KML document to a file
        kml_buffer = BytesIO()
        kml_buffer.write(k.to_string(prettyprint=True).encode('utf-8'))
        kml_buffer.seek(0)
        st.write("Done.")

        st.download_button(
            label="Download RT/RW Net (Telkomsel) KML",
            data=kml_buffer,
            file_name="Telkomsel_RT_RW_Net.kml",
            mime="application/vnd.google-earth.kml+xml"
        )

    elif export_option == "RT/RW (Others)":
        # Initialize KML object
        st.write("Exporting RT/RW Net (Others) file...")
        k = kml.KML()
        
        # Create a KML document
        document = kml.Document()
        k.append(document)
        
        # Create a KML placemark for each point
        for _, row in st.session_state.others_rt_rw_net_df.iterrows():
            placemark = kml.Placemark()
            placemark.geometry = Point(row["longitude"], row["latitude"])
            placemark.name = f"{row['provider']} (Lat: {row['latitude']}, Lon: {row['longitude']})"
            
            # Add attributes as KML extended data
            extended_data = kml.ExtendedData()
            extended_data.elements.append(kml.Data(name="desa", value=str(row["desa"])))
            extended_data.elements.append(kml.Data(name="kecamatan", value=str(row["kecamatan"])))
            extended_data.elements.append(kml.Data(name="kabupaten", value=str(row["kabupaten"])))
            extended_data.elements.append(kml.Data(name="provider", value=str(row["provider"])))
            extended_data.elements.append(kml.Data(name="latency", value=str(row["latency"])))
            placemark.extended_data = extended_data
            
            document.append(placemark)
        
        # Save the KML document to a file
        kml_buffer = BytesIO()
        kml_buffer.write(k.to_string(prettyprint=True).encode('utf-8'))
        kml_buffer.seek(0)
        st.write("Done.")

        st.download_button(
            label="Download RT/RW Net (Others) KML",
            data=kml_buffer,
            file_name="Others_RT_RW_Net.kml",
            mime="application/vnd.google-earth.kml+xml"
        )

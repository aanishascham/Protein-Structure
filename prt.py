import streamlit as st
from stmol import showmol
import py3Dmol
import requests
from Bio.PDB import PDBParser
from Bio import PDB
from difflib import SequenceMatcher
from io import StringIO
import io
import plotly.graph_objects as go

st.set_page_config(
    page_title='ProteinVista', 
    page_icon="🧬", 
    layout='wide', 
    initial_sidebar_state='auto'
)

# Your Streamlit app code here
st.sidebar.title('Structure Spin & PDB Sync')
st.sidebar.write('''With our innovative protein visualization app,
Input your structure, see it vividly unwrap,
Interact with the model, explore and spin,
Download the PDB file, where discoveries begin!''')
st.markdown(
    r"""
    <style>
    .stDeployButton {
        visibility: hidden;
    }
    </style>
    """,
    unsafe_allow_html=True
)

st.markdown(
    """
    <style>
    .stMultiSelect.clearable div[role="option"]:first-child,
    .stMultiSelect.clearable div[role="option"]:last-child {
        display: none !important;
    }
    </style>
    """,
    unsafe_allow_html=True
)

# stmol
def render_mol(pdb):
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb, 'pdb')
    pdbview.setStyle({'cartoon': {'color': 'spectrum'}})
    pdbview.setBackgroundColor('white')
    pdbview.zoomTo()
    pdbview.zoom(2, 800)
    pdbview.spin(True)
    showmol(pdbview, height=500, width=800)

# Protein sequence input
DEFAULT_SEQ = "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"
txt = st.sidebar.text_area('Input sequence', DEFAULT_SEQ, height=275)



# ESMfold
def update(sequence=txt):
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
    }
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence, verify=False)

    name = sequence[:3] + sequence[-3:]
    pdb_string = response.content.decode('utf-8')

    with open('predicted.pdb', 'w') as f:
        f.write(pdb_string)

    parser = PDBParser()
    struct = parser.get_structure('predicted', 'predicted.pdb')
    atoms_list = list(struct.get_atoms())
    b_value = round(sum(atom.get_bfactor() for atom in atoms_list) / len(atoms_list), 4)


    # Display protein structure
    st.subheader('Visualization of predicted protein structure')
    render_mol(pdb_string)

    # plDDT value is stored in the B-factor field
    st.subheader('plDDT')
    st.write('plDDT is a per-residue estimate of the confidence in prediction')
    st.info(f'plDDT: {b_value}')

    st.download_button(
        label="Download PDB",
        data=pdb_string,
        file_name='predicted.pdb',
        mime='text/plain',
    )

predict = st.sidebar.button('Predict', on_click=update)


if not predict:
    st.warning('Please input protein sequence data.')
st.write("---") 


st.sidebar.write('''Compare PDB Files''')
def compare_pdb_files(file1, file2):
    content1 = file1.read()
    content2 = file2.read()

    lines1 = content1.decode("utf-8").split('\n')
    lines2 = content2.decode("utf-8").split('\n')

    # Calculate the similarity using the Ratcliff/Obershelp algorithm
    matcher = SequenceMatcher(None, lines1, lines2)
    similarity_ratio = matcher.ratio()
    similarity_percent = similarity_ratio * 100

    # Generate a list of differences
    differences = list(matcher.get_opcodes())

    return similarity_percent, differences

def main():
    st.title("PDB File Comparator")

    st.sidebar.title("Upload PDB Files")
    uploaded_file1 = st.sidebar.file_uploader("Upload PDB File 1", type="pdb")
    uploaded_file2 = st.sidebar.file_uploader("Upload PDB File 2", type="pdb")

    if uploaded_file1 and uploaded_file2:
        st.subheader("Uploaded PDB Files")
        st.write("PDB File 1:", uploaded_file1.name)
        st.write("PDB File 2:", uploaded_file2.name)

        similarity_percent, differences = compare_pdb_files(uploaded_file1, uploaded_file2)
        
        st.subheader("Similarity Percentage")
        st.write(f"The files match {similarity_percent:.2f}%")

        if differences:
            st.subheader("Differences Found")
            st.write("Index", "Operation", "File 1", "File 2")
            lines1 = uploaded_file1.getvalue().decode("utf-8").split('\n')
            lines2 = uploaded_file2.getvalue().decode("utf-8").split('\n')
            for diff in differences:
                operation, start1, end1, start2, end2 = diff
                if operation != 'equal':
                    st.write(start1, operation, lines1[start1:end1], lines2[start2:end2])

        else:
            st.subheader("No Differences Found")

if __name__ == "__main__":
    main()

st.write("---") 

# Sidebar title and file uploader
st.sidebar.title('Upload PDB File that you want to read')
uploaded_file = st.sidebar.file_uploader("Choose a PDB file", type=['pdb'])

st.markdown(
    r"""
    <style>
    .stDeployButton {
        visibility: hidden;
    }
    </style>
    """,
    unsafe_allow_html=True
)

st.markdown(
    """
    <style>
    .stMultiSelect.clearable div[role="option"]:first-child,
    .stMultiSelect.clearable div[role="option"]:last-child {
        display: none !important;
    }
    </style>
    """,
    unsafe_allow_html=True
)

if uploaded_file is not None:
    # Read the contents of the uploaded PDB file
    pdb_string = uploaded_file.read().decode('utf-8')

    # Display the raw text of the PDB file
    st.subheader('PDB File Contents')
    st.code(pdb_string)
st.write("---") 

def get_ligand_ids(pdb_content):
    pdb_content_string = io.StringIO(pdb_content.decode('utf-8'))
    parser = PDB.PDBParser()
    structure = parser.get_structure('protein', pdb_content_string)

    ligand_ids = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                hetero_flag = residue.get_id()[0]
                if hetero_flag[0] != ' ':
                    ligand_ids.add(residue.id[1])
    return sorted(list(ligand_ids))

def visualize_binding_sites(pdb_content, ligand_id):
    pdb_content_string = io.StringIO(pdb_content.decode('utf-8'))
    parser = PDB.PDBParser()
    structure = parser.get_structure('protein', pdb_content_string)

    ligand_coordinates = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[1] == ligand_id:
                    for atom in residue:
                        ligand_coordinates.append(atom.get_coord())

    if not ligand_coordinates:
        st.error("Ligand ID not found in the protein structure.")
        return

    protein_coordinates = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    protein_coordinates.append(atom.get_coord())

    protein_x, protein_y, protein_z = zip(*protein_coordinates)
    ligand_x, ligand_y, ligand_z = zip(*ligand_coordinates)

    fig = go.Figure()

    # Add protein trace
    fig.add_trace(go.Scatter3d(
        x=protein_x, y=protein_y, z=protein_z,
        mode='markers', marker=dict(color='blue', size=3),
        name='Protein'
    ))

    # Add ligand trace
    fig.add_trace(go.Scatter3d(
        x=ligand_x, y=ligand_y, z=ligand_z,
        mode='markers', marker=dict(color='red', size=5),
        name='Ligand'
    ))

    fig.update_layout(
        scene=dict(
            xaxis_title='X Axis',
            yaxis_title='Y Axis',
            zaxis_title='Z Axis',
        ),
        title='Protein-Ligand Binding Sites Visualization',
        margin=dict(l=0, r=0, t=40, b=0),
    )

    st.plotly_chart(fig)

def main():
    st.title("Protein-Ligand Binding Sites Visualization")

    pdb_file = st.file_uploader("Upload PDB file", type=["pdb"])

    if pdb_file is not None:
        st.write("Uploaded PDB file:", pdb_file.name)
        
        ligand_ids = get_ligand_ids(pdb_file.read())
        if ligand_ids:
            selected_ligand_id = st.selectbox("Select Ligand ID:", ligand_ids)

            if st.button("Visualize"):
                pdb_content = pdb_file.getvalue()
                visualize_binding_sites(pdb_content, selected_ligand_id)
        else:
            st.warning("No ligand IDs found in the PDB file.")

if __name__ == "__main__":
    main()

const DemographicModal = ({ show, setShow }) => {

  return (
    <div>
      <Modal open={show} onClose={setShow}>
        <div>
          <h2>{demographic.name}</h2>
          <p>{demographic.age}</p>
        </div>
      </Modal>
    </div>
  );
}

export default DemographicModal;
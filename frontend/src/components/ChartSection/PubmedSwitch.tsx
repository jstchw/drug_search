import React from "react";
import { Form } from "react-bootstrap";

interface PubmedSwitchProps {
  handlePubmedSwitch: () => void;
}

const PubmedSwitch: React.FC<PubmedSwitchProps> = ({ handlePubmedSwitch }) => {

  return (
    <Form className="d-flex justify-content-center">
      <Form.Check
        type="switch"
        id="custom-switch"
        label="Include Pubmed data"
        onChange={handlePubmedSwitch}
      />
    </Form>
  );
};

export default PubmedSwitch;

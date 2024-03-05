import React from 'react';
import { Form } from 'react-bootstrap';

interface ToggleSwitchProps {
  handleToggleSwitch: () => void;
  children: React.ReactNode;
}

const ToggleSwitch: React.FC<ToggleSwitchProps> = ({ handleToggleSwitch, children }) => {
  return (
    <Form className="d-flex justify-content-center">
      <Form.Check type="switch" id="custom-switch" onChange={handleToggleSwitch} />
      <Form.Check.Label>{children}</Form.Check.Label>
    </Form>
  );
};

export default ToggleSwitch;

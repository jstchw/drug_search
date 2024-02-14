import React from "react";
import { Toast, ToastContainer, Col } from "react-bootstrap";
import { useUrlParams } from "../../hooks/useUrlParams";
import { URLParams } from "../../types";
import { mapParamArrayToLabels } from "../../utils/utils";

const InfoCard = () => {
  const { params, paramError } = useUrlParams();
  const show = !paramError;

  const constraints = mapParamArrayToLabels(params);

  return (
    <ToastContainer
      className={"patient-card m-3"}
      position={"bottom-end"}
      containerPosition={"fixed"}
    >
      <Toast show={show} delay={3000} style={{ width: "fit-content" }}>
        <Toast.Header className={"justify-content-end"} closeButton={false}>
          <span>Info Card</span>
        </Toast.Header>
        <Toast.Body>
          {Object.keys(constraints).map((key) => {
            const value = constraints[key as keyof URLParams];
            return (
              <React.Fragment key={key}>
                {value && (
                  <Col
                    className={
                      "fs-5 d-flex justify-content-start align-items-center"
                    }
                  >
                    <span className={"text-secondary"}>{key}</span>
                    <div className={"vr mx-2"} />
                    <span>{value}</span>
                  </Col>
                )}
              </React.Fragment>
            );
          })}
        </Toast.Body>
      </Toast>
    </ToastContainer>
  );
};

export default InfoCard;

import React, { ReactNode } from "react";
import { Button } from "react-bootstrap";
import "./ReadMore.css";

type ReadMoreProps = {
  children: ReactNode;
};

export const ReadMore: React.FC<ReadMoreProps> = ({ children }) => {
  const [isExpanded, setIsExpanded] = React.useState(true);
  const childrenArray = React.Children.toArray(children);
  const isTooLong = childrenArray.length > 20;

  const buttonText = isExpanded ? "show more" : "show less";
  return (
    <React.Fragment>
      <p>{isExpanded ? childrenArray.slice(0, 20) : childrenArray}</p>
      <div className={"d-flex justify-content-center"}>
        {isTooLong && (
          <Button
            size={"sm"}
            className={"expand-button"}
            onClick={() => setIsExpanded(!isExpanded)}
          >
            {buttonText}
          </Button>
        )}
      </div>
    </React.Fragment>
  );
};

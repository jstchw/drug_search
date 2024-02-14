import React from "react";
import { Popover, OverlayTrigger, Badge } from "react-bootstrap";
import { useUrlParams } from "../../hooks/useUrlParams";
import { drugGroupConfig } from "../../constants";

interface DrugGroupsProps {
  drugGroups: string[];
}

const DrugGroups: React.FC<DrugGroupsProps> = ({ drugGroups }) => {
  const { params } = useUrlParams();
  const processDrugGroups = (groups: string[]) => {
    return groups.map((group, index) => {
      const groupData = drugGroupConfig[group];

      if (!groupData) {
        throw new Error(`No data found for group ${group}`);
      }

      const { IconComponent, label, variant, description } = groupData;
      const groupLabel = group.substring(0, 3).toUpperCase();

      const popover = (
        <Popover id={`popover-${index}`}>
          <Popover.Header style={{ backgroundColor: "transparent" }}>
            <div className={"d-flex align-items-center fs-5"}>
              {IconComponent && <IconComponent weight={"light"} />}
              <div className={"vr mx-2"} />
              {label}
            </div>
          </Popover.Header>
          <Popover.Body>{description}</Popover.Body>
        </Popover>
      );

      return (
        <OverlayTrigger
          key={index}
          trigger={["hover", "focus"]}
          placement="bottom"
          overlay={popover}
        >
          <Badge style={{ cursor: "default" }} className="mx-1" bg={variant}>
            {groupLabel}
          </Badge>
        </OverlayTrigger>
      );
    });
  };

  return (
    <>
      {params.searchBy === "generic_name" ? (
        <h3>{processDrugGroups(drugGroups)}</h3>
      ) : (
        <h5>{processDrugGroups(drugGroups)}</h5>
      )}
    </>
  );
};

export default DrugGroups;

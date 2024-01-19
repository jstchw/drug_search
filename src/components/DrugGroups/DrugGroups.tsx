import React from 'react';
import { Popover, OverlayTrigger, Badge } from 'react-bootstrap';
import { IconProps, Prescription, MagnifyingGlass, PawPrint, Trash, Carrot, Flask, X, SmileyNervous } from "@phosphor-icons/react";

interface DrugGroupConfig {
    [key: string]: {
        label: string;
        variant: string;
        description: string;
        IconComponent: React.ForwardRefExoticComponent<IconProps>
    };
}

const groupConfig: DrugGroupConfig = {
    approved: {
        label: 'Approved',
        variant: 'success',
        description: 'Has been approved in at least one jurisdiction, at some point in time.',
        IconComponent: Prescription
    },
    investigational: {
        label: 'Investigational',
        variant: 'warning',
        description: 'Undergoing evaluation in the drug approval process in at least one jurisdiction.',
        IconComponent: MagnifyingGlass
        },
    illicit: {
        label: 'Illicit',
        variant: 'danger',
        description: 'Deemed illegal in at least one jurisdiction, at some point in time.',
        IconComponent: X
    },
    vet_approved: {
        label: 'Vet Approved',
        variant: 'info',
        description: 'Approved for veterinary use in at least one jurisdiction, at some point in time.',
        IconComponent: PawPrint
    },
    withdrawn: {
        label: 'Withdrawn',
        variant: 'secondary',
        description: 'Has been withdrawn from the market in at least one jurisdiction, at some point in time.',
        IconComponent: Trash
    },
    nutraceutical: {
        label: 'Nutraceutical',
        variant: 'light',
        description: 'Pharmaceutical-grade nutrient with potential health benefits.',
        IconComponent: Carrot
    },
    experimental: {
        label: 'Experimental',
        variant: 'dark',
        description: 'Shown to bind specific proteins in experimental settings.',
        IconComponent: Flask
    },
    side_effect: {
        label: 'Side Effect',
        variant: 'primary',
        description: 'Has been reported as a side effect.',
        IconComponent: SmileyNervous
    }
};

interface DrugGroupsProps {
    drugGroups: string[];
}

const DrugGroups: React.FC<DrugGroupsProps> = ({ drugGroups }) => {
    const processDrugGroups = (groups: string[]) => {
        return groups.map((group, index) => {
            const groupData = groupConfig[group];

            if (!groupData) {
                throw new Error(`No data found for group ${group}`);
            }

            const { IconComponent, label, variant, description } = groupData;
            const groupLabel = group.substring(0, 3).toUpperCase();

            const popover = (
                <Popover id={`popover-${index}`}>
                    <Popover.Header
                        style={{backgroundColor: 'transparent'}}
                    >
                        <div className={'d-flex align-items-center fs-5'}>
                            {IconComponent && <IconComponent weight={'light'}/>}
                            <div className={'vr mx-2'}/>
                            {label}
                        </div>
                    </Popover.Header>
                    <Popover.Body>{description}</Popover.Body>
                </Popover>
            );

            return (
                <OverlayTrigger key={index} trigger={['hover', 'focus']} placement="bottom" overlay={popover}>
                    <Badge style={{cursor: "default"}} className="mx-1" bg={variant}>{groupLabel}</Badge>
                </OverlayTrigger>
            );
        });
    };

    return (
        <h3>{processDrugGroups(drugGroups)}</h3>
    );
};

export default DrugGroups;

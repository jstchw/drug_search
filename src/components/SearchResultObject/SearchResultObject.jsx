import React, {useEffect} from "react";
import {Col, Container, Row, Badge, Placeholder, Popover, OverlayTrigger, Alert, Button} from "react-bootstrap";
import {Book as BookIcon} from "react-bootstrap-icons";
import './SearchResultObject.css'
import DrugDescription from "../DrugDescription/DrugDescription";
import DrugAccordion from "../DrugAccordion/DrugAccordion";
import ChartDisplayObject from "../ChartDisplayObject/ChartDisplayObject";
import PatientCard from "../PatientCard/PatientCard";
import {searchTypes} from "../OptionModal/OptionModal";
import useDrugInfo from "../../hooks/useDrugInfo";
import DemographicModal from "../DemographicModal/DemographicModal";
import { isMobile } from 'react-device-detect'

const getUniqueObjects = (array) => {
    return array.filter((obj, index, self) =>
        index === self.findIndex((otherObj) => obj.drug_name === otherObj.drug_name)
    );
}


const SearchResultObject = (props) => {
    const { termCountDict, totalCount } = props.searchResults.result
    const eventsOverTime = props.eventsOverTime.result.results
    const searchTypeRef = React.useRef(props.searchOptions.searchBy)

    const drugInfo = useDrugInfo(props.searchTerm, searchTypeRef.current)

    const [showDemographicModal, setShowDemographicModal] = React.useState(false)
    const [selectedWord, setSelectedWord] = React.useState(null)

    useEffect(() => {
        if(!showDemographicModal) {
            setSelectedWord(null);
        }
    }, [showDemographicModal])

    const uniqueDrugInfo = getUniqueObjects(drugInfo)
    let groupDescription

    const processDrugGroups = (groups) => {
        return groups.map((group, index) => {
            let variant
            switch (group) {
                case 'approved':
                    variant = 'success'
                    groupDescription = 'Has been approved in at least one jurisdiction, at some point in time.'
                    break
                case 'investigational':
                    variant = 'warning'
                    groupDescription = 'Undergoing evaluation in the drug approval process in at least one jurisdiction.'
                    break
                case 'illicit':
                    variant = 'danger'
                    groupDescription = 'Deemed illegal in at least one jurisdiction, at some point in time.'
                    break
                case 'vet_approved':
                    variant = 'info'
                    groupDescription = 'Approved for veterinary use in at least one jurisdiction, at some point in time.'
                    break
                case 'withdrawn':
                    variant = 'secondary'
                    groupDescription = 'Has been withdrawn from the market in at least one jurisdiction, at some point in time.'
                    break
                case 'nutraceutical':
                    variant = 'light'
                    groupDescription = 'Pharmaceutical-grade nutrient with potential health benefits.'
                    break
                case 'experimental':
                    variant = 'dark'
                    groupDescription = 'Shown to bind specific proteins in experimental settings.'
                    break
                default:
                    variant = 'secondary'
            }
            group = group.substring(0, 3).toUpperCase()

            const popover = (
                <Popover id={`popover-${index}`}>
                    <Popover.Body>
                        <div className={'fs-6 text-center'}>{groupDescription}</div>
                    </Popover.Body>
                </Popover>
            )

            return (
                <OverlayTrigger key={index} trigger={['hover', 'focus']} placement={'bottom'} overlay={popover}>
                    <Badge className={'mx-1 no-pointer-badge'} key={index} bg={variant}>{group}</Badge>
                </OverlayTrigger>
            )
        })

    }

    const groupedByProduct = drugInfo.reduce((acc, term) => {
        if (!acc[term.product]) {
            acc[term.product] = [];
        }
        acc[term.product].push(term);
        return acc;
    }, {});

    const DrugInfoObject = ({drugInfo}) => {
        return (
            drugInfo.length > 0 ? (
                Object.entries(groupedByProduct).map(([product, terms], index) => (
                    <Col key={index}>
                        {/* Searching by active substance*/}
                        {searchTypeRef.current === searchTypes[0].value &&
                            <>
                                <h1>{terms[0].drug_name}</h1>
                                {terms[0].groups && <h3>{processDrugGroups(terms[0].groups)}</h3>}
                                {/* Demographic modal for every drug */}
                                <Button onClick={() => {
                                    setSelectedWord(terms[0].drug_name)
                                    setShowDemographicModal(true)
                                }}>
                                    <BookIcon size={'22'}></BookIcon>
                                </Button>
                            </>
                        }
                        {/* Searching by brand names*/}
                        {searchTypeRef.current === searchTypes[1].value &&
                            <>
                                <h1>{product}</h1>
                                <Badge className={'mb-3'}>Active substance</Badge>
                                {terms.map((term, termIndex) => (
                                    <React.Fragment key={termIndex}>
                                        <h5>
                                            {term.drug_name}{processDrugGroups(term.groups)}
                                            <Button className={'mx-1'} size={'sm'} onClick={() => {
                                                setSelectedWord(term.drug_name)
                                                setShowDemographicModal(true)
                                            }}>
                                                <BookIcon/>
                                            </Button>
                                        </h5>
                                    </React.Fragment>
                                ))}
                            </>
                        }
                        {/* Searching by adverse effect*/}
                        {searchTypeRef.current === searchTypes[2].value &&
                            <>
                                <h1>{terms.map(term => term.ADE).join(' & ')}</h1>
                            </>
                        }
                    </Col>
                ))
            ) : (
                <>
                    <Placeholder as={'h1'} animation="glow">
                        <Placeholder xs={6} />
                    </Placeholder>
                    <Placeholder as={'h3'} animation="glow">
                        <Placeholder xs={3} />
                    </Placeholder>
                </>
            )
        )
    }

    return (
        <div>
            <Container className={'my-4'}>
                {!isMobile ? <Row><DrugInfoObject drugInfo={drugInfo}/></Row> : <DrugInfoObject drugInfo={drugInfo}/>}
            </Container>

            <Container className={'main-drug-container'}>
                <Alert variant={"secondary"}>
                    The provided information is indicative and shouldn't be used for inference.
                </Alert>
                { !isMobile  &&
                    <Row>
                        {(uniqueDrugInfo.length > 0 && uniqueDrugInfo
                                .every((term) => term.ADE === null && term.full_info === true)) &&
                            (uniqueDrugInfo.map((term, index) => (
                                uniqueDrugInfo.length > 1 ? (
                                    <Col key={index}>
                                        <DrugDescription drugInfo={term} showAdditionalSearch={props.showAdditionalSearch}/>
                                        <DrugAccordion drugInfo={term} totalCount={totalCount}/>
                                    </Col>
                                ) : (
                                    <React.Fragment key={index}>
                                        <Col>
                                            <DrugDescription drugInfo={term} showAdditionalSearch={props.showAdditionalSearch}/>
                                        </Col>
                                        <Col>
                                            <DrugAccordion drugInfo={term} totalCount={totalCount}/>
                                        </Col>
                                    </React.Fragment>
                                )
                            ))
                        )}
                    </Row>
                }
                { isMobile &&
                    <React.Fragment>
                        <Col>
                            {(uniqueDrugInfo.length > 0 && uniqueDrugInfo
                                .every((term) => term.ADE === null && term.full_info === true)) &&
                            (uniqueDrugInfo.map((term, index) => (
                                    <React.Fragment key={index}>
                                        <div className={'mb-4'}>
                                            <DrugDescription drugInfo={term}
                                                             showAdditionalSearch={props.showAdditionalSearch}/>
                                            <DrugAccordion drugInfo={term} totalCount={totalCount}/>
                                        </div>
                                    </React.Fragment>
                                ))
                            )}
                        </Col>
                    </React.Fragment>
                }

            </Container>

            <Col>
                <ChartDisplayObject eventsOverTime={eventsOverTime} termCountDict={termCountDict}
                                    totalCount={totalCount} searchResults={termCountDict}
                                    searchOptions={props.searchOptions}
                />
            </Col>

            <DemographicModal
                show={showDemographicModal}
                handleClose={() => setShowDemographicModal(false)}
                selectedWord={selectedWord}
                isMobile={isMobile}
            />
            { !isMobile && <PatientCard searchOptions={props.searchOptions} totalADE={totalCount} /> }
        </div>
    );
}
export default SearchResultObject;